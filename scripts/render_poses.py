#!/usr/bin/env python3
"""
scripts/render_poses.py — 도킹 포즈를 PNG 이미지로 렌더링한다.

3Dmol.js 를 헤드리스 Chromium 에서 띄워 스크린샷을 찍는다. 실험마다 포즈를
사진으로 남겨야 숫자만 보고 놓치는 것을 잡을 수 있다. EPAS1 에서 EF1% 32.68 이
좋아 보였지만 상위 히트 그림을 보고서야 전부 같은 화학형임을 알았다.

출력
  poses/<화합물>.png        수용체 카툰 + 포켓 표면 + 리간드 스틱
  poses/_overlay.png        공결정 리간드 자리에 벤치마크들을 겹쳐 본 그림
"""
from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
THREEDMOL = ROOT / "node_modules" / "3dmol" / "build" / "3Dmol-min.js"

PAGE = """<!doctype html>
<html><head><meta charset="utf-8"><style>
  html,body{{margin:0;padding:0;background:#0d1117;overflow:hidden}}
  #v{{width:{w}px;height:{h}px;position:relative}}
  #cap{{position:absolute;left:14px;bottom:12px;z-index:10;
       font:13px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#e6edf3;
       background:rgba(13,17,23,.82);padding:8px 12px;border-radius:6px;
       border:1px solid #30363d;white-space:pre}}
</style><script>{lib}</script></head>
<body><div id="v"></div><div id="cap">{caption}</div>
<script>
const rec = {rec};
const lig = {lig};
const extra = {extra};
const v = $3Dmol.createViewer(document.getElementById('v'),
          {{backgroundColor:'#0d1117'}});

// 수용체 먼저, 그다음 리간드 — 모델 객체를 잡아야 근접 선택이 된다
const rm = v.addModel(rec,'pdb');
const lm = v.addModel(lig,'sdf');

v.setStyle({{model:rm}}, {{cartoon:{{color:'#56607a', opacity:0.72}}}});
// 리간드 5 A 이내 잔기만 스틱으로 — 포켓을 이루는 잔기들
v.addStyle({{model:rm, byres:true, within:{{distance:5.0, sel:{{model:lm}}}}}},
           {{stick:{{radius:0.11, colorscheme:'grayCarbon'}}}});
v.setStyle({{model:lm}}, {{stick:{{radius:0.26, colorscheme:'cyanCarbon'}}}});

extra.forEach((e) => {{
  const em = v.addModel(e.sdf,'sdf');
  v.setStyle({{model:em}}, {{stick:{{radius:0.18, color:e.color}}}});
}});

v.zoomTo({{model:lm}});
v.zoom(0.42);   // 뒤로 빼서 포켓 맥락이 보이게
v.render();
window.__ready = true;
</script></body></html>
"""


def mol_block(mol: Chem.Mol) -> str:
    return Chem.MolToMolBlock(mol)


def load_poses(sdf: Path) -> dict:
    out = {}
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False):
        if m is None:
            continue
        name = m.GetProp("_Name").rsplit("_pose", 1)[0] if m.HasProp("_Name") else "?"
        if name in out:
            continue
        out[name] = {
            "mol": m,
            "group": m.GetProp("group") if m.HasProp("group") else "unknown",
            "score": float(m.GetProp("vina_score")) if m.HasProp("vina_score") else None,
        }
    return out


def render(page_html: str, png: Path, w: int, h: int) -> bool:
    from playwright.sync_api import sync_playwright

    tmp = png.with_suffix(".html")
    tmp.write_text(page_html)
    try:
        with sync_playwright() as p:
            b = p.chromium.launch(args=["--no-sandbox", "--disable-gpu", "--use-gl=swiftshader"])
            pg = b.new_page(viewport={"width": w, "height": h}, device_scale_factor=2)
            pg.goto(f"file://{tmp.resolve()}")
            pg.wait_for_function("window.__ready === true", timeout=60000)
            pg.wait_for_timeout(2500)  # 표면 렌더 완료 대기
            pg.locator("#v").screenshot(path=str(png))
            b.close()
        return True
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expdir", required=True, help="exp/<GENE> 디렉터리")
    ap.add_argument("--arm", default="treatment")
    ap.add_argument("--top", type=int, default=6, help="상위 히트 몇 개를 렌더링할지")
    ap.add_argument("--width", type=int, default=900)
    ap.add_argument("--height", type=int, default=700)
    args = ap.parse_args()

    expdir = Path(args.expdir)
    armdir = expdir / args.arm
    outdir = expdir / "poses"
    outdir.mkdir(parents=True, exist_ok=True)

    lib = THREEDMOL.read_text()
    rec = (armdir / "receptor.pdb").read_text()
    meta = json.loads((armdir / "target_meta.json").read_text())
    center = meta["pocket"]["center"]
    pdb_id = meta.get("pdb_id")

    poses = load_poses(armdir / "poses.sdf")
    bench = [(k, v) for k, v in poses.items() if v["group"] == "benchmark"]
    hits = sorted(
        ((k, v) for k, v in poses.items() if v["group"] == "top_hit" and v["score"]),
        key=lambda kv: kv[1]["score"],
    )[: args.top]

    made = []
    for name, r in bench + hits:
        cap = (
            f"{name}\n"
            f"Vina {r['score']:.2f} kcal/mol   {r['group']}\n"
            f"{pdb_id}  pocket {center}"
        )
        html = PAGE.format(
            w=args.width,
            h=args.height,
            lib=lib,
            caption=cap,
            rec=json.dumps(rec),
            lig=json.dumps(mol_block(r["mol"])),
            extra="[]",
            center=json.dumps(center),
        )
        png = outdir / f"{name}.png"
        print(f"  렌더링 {name} ... ", end="", flush=True)
        ok = render(html, png, args.width, args.height)
        print("OK" if ok and png.exists() else "실패")
        if png.exists():
            made.append({"name": name, "group": r["group"], "score": r["score"],
                         "png": str(png.relative_to(expdir))})

    # 벤치마크 3개를 한 그림에 겹쳐 — 같은 포켓에 어떻게 앉는지
    if len(bench) >= 2:
        colors = ["#f778ba", "#d29922", "#3fb950"]
        first = bench[0][1]["mol"]
        extra = [
            {"sdf": mol_block(v["mol"]), "color": colors[i % len(colors)]}
            for i, (k, v) in enumerate(bench[1:])
        ]
        cap = "벤치마크 중첩\n" + "   ".join(
            f"{k} {v['score']:.2f}" for k, v in bench
        )
        html = PAGE.format(
            w=args.width, h=args.height, lib=lib, caption=cap,
            rec=json.dumps(rec), lig=json.dumps(mol_block(first)),
            extra=json.dumps(extra), center=json.dumps(center),
        )
        png = outdir / "_benchmark_overlay.png"
        print("  렌더링 벤치마크 중첩 ... ", end="", flush=True)
        ok = render(html, png, args.width, args.height)
        print("OK" if ok and png.exists() else "실패")
        if png.exists():
            made.append({"name": "benchmark_overlay", "group": "overlay",
                         "score": None, "png": str(png.relative_to(expdir))})

    (expdir / "poses" / "index.json").write_text(
        json.dumps(made, indent=2, ensure_ascii=False)
    )
    print(f"\n[+] 이미지 {len(made)}장 → {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
