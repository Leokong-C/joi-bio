#!/usr/bin/env python3
"""
scripts/render_m252.py — M252 측쇄 뒤집힘을 before/after 로 그린다.

무엇을 보여주는가
  Ren et al. 2022 (Mol Pharmacol 102:240) 은 belzutifan 이 HIF-2α 에 결합하면
  MET252 측쇄가 포켓 밖으로 뒤집힌다고 보고했다. 우리 리지드 도킹(6CZW)에서
  M252 는 belzutifan 과 1.53~2.05 Å 거리에 있었다 — 탄소-탄소 반데르발스 접촉이
  3.4~4.0 Å 이므로 충돌이다. 측쇄를 풀어 다시 도킹하자 M252 가 실제로 회전했다.

    MET252  측쇄 RMSD 2.12 Å, 말단 CE 3.96 Å 이동, 접촉 1.53 → 3.21 Å
    PHE254  0.15 Å — 자유를 줬으나 움직이지 않음
    HIS293  0.37 Å — 움직이지 않음

  셋 중 하나만 움직였다는 점이 핵심이다. '자유도를 늘려 점수가 좋아진' 잡음이
  아니라 특정 잔기의 특정 회전이라는 뜻이다. 그림이 그걸 보여줘야 한다.

입력
  exp/<gene>/treatment/receptor.pdb          원래 구조 (6CZW 체인 A)
  exp/<gene>/flex/poses/<chembl>.pdbqt       유연 도킹 출력 (리간드 + 움직인 측쇄)
"""
from __future__ import annotations

import argparse
import base64
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THREEDMOL = ROOT / "node_modules" / "3dmol" / "build" / "3Dmol-min.js"

ELEM = {"C": "C", "N": "N", "O": "O", "S": "S", "F": "F", "H": "H",
        "A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H"}


def parse_flex_pose(path: Path, model: int = 1) -> tuple[list, dict]:
    """유연 도킹 PDBQT 의 MODEL n 에서 리간드 원자와 유연 잔기를 뽑는다.

    구조: MODEL n → (리간드 ROOT/BRANCH) → BEGIN_RES <잔기> ... END_RES → ENDMDL
    """
    lig: list[tuple[str, tuple[float, float, float], str]] = []
    flex: dict[str, list] = {}
    cur = None
    inside = False
    for line in path.read_text().splitlines():
        if line.startswith(f"MODEL {model}"):
            inside = True
            continue
        if inside and line.startswith("ENDMDL"):
            break
        if not inside:
            continue
        if line.startswith("BEGIN_RES"):
            p = line.split()
            cur = (p[1], p[2], int(p[3]))  # (resname, chain, resid)
            flex[f"{p[1]}{p[3]}"] = []
            continue
        if line.startswith("END_RES"):
            cur = None
            continue
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            # PDBQT 는 자체 원자타입을 쓴다. 뒤쪽 필드에서 원소를 복원한다.
            at = line[77:79].strip() or line[76:78].strip()
            el = ELEM.get(at.upper(), name[0])
            (flex[f"{cur[0]}{cur[2]}"] if cur else lig).append((name, xyz, el))
    return lig, flex


def residue_atoms(pdb: Path, resname: str, resid: int, chain: str = "A") -> list:
    out = []
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if (line[17:20].strip() == resname and line[22:26].strip() == str(resid)
                and line[21] == chain):
            out.append((line[12:16].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                        (line[76:78].strip() or line[12:16].strip()[0])))
    return out


def to_pdb(atoms: list, resname: str, resid: int, het: bool = False) -> str:
    rec = "HETATM" if het else "ATOM  "
    lines = []
    for i, (name, (x, y, z), el) in enumerate(atoms, 1):
        nm = f" {name:<3s}" if len(name) < 4 else name
        lines.append(
            f"{rec}{i:5d} {nm}{resname:>4s} A{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}"
        )
    lines.append("END")
    return "\n".join(lines)


def centroid(atoms: list) -> tuple[float, float, float]:
    n = len(atoms)
    return tuple(sum(a[1][i] for a in atoms) / n for i in range(3))  # type: ignore


def closest_pair(a: list, b: list) -> tuple[float, tuple, tuple]:
    best = (1e9, None, None)
    for _, p, _ in a:
        for _, q, _ in b:
            d = math.dist(p, q)
            if d < best[0]:
                best = (d, p, q)
    return best  # type: ignore


PAGE = """<!doctype html>
<html><head><meta charset="utf-8"><style>
  html,body{margin:0;padding:0;background:#0d1117;overflow:hidden;
            font:13px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#e6edf3}
  .wrap{display:flex;width:__TOTW__px;height:__H__px}
  .pane{position:relative;width:__W__px;height:__H__px}
  .pane + .pane{border-left:1px solid #30363d}
  .ttl{position:absolute;left:14px;top:12px;z-index:10;font-size:15px;font-weight:600}
  .sub{position:absolute;left:14px;top:34px;z-index:10;color:#8b949e}
  .cap{position:absolute;left:14px;bottom:12px;z-index:10;white-space:pre;
       background:rgba(13,17,23,.82);padding:8px 12px;border-radius:6px;
       border:1px solid #30363d}
  .key{position:absolute;right:14px;top:12px;z-index:10;white-space:pre;text-align:right}
</style><script>__LIB__</script></head>
<body>
<div class="wrap">
  <div class="pane" id="p0">
    <div class="ttl">__T0__</div><div class="sub">__S0__</div>
    <div class="key" id="k0"></div><div class="cap">__C0__</div>
  </div>
  <div class="pane" id="p1">
    <div class="ttl">__T1__</div><div class="sub">__S1__</div>
    <div class="key" id="k1"></div><div class="cap">__C1__</div>
  </div>
</div>
<script>
const D = __DATA__;

function build(paneId, keyId, spec) {
  const v = $3Dmol.createViewer(document.getElementById(paneId),
            {backgroundColor:'#0d1117'});
  const rm = v.addModel(D.receptor, 'pdb');
  v.setStyle({model:rm}, {cartoon:{color:'#657090', opacity:0.40}});

  const lm = v.addModel(D.ligand, 'pdb');
  v.setStyle({model:lm}, {stick:{radius:0.24, colorscheme:'cyanCarbon'}});

  // 포켓 맥락 — 리간드 5 A 이내 잔기를 가는 스틱으로. 선(line)은 스크린샷에서
  // 거의 보이지 않아 포켓이 어디인지 읽히지 않았다.
  v.addStyle({model:rm, byres:true, within:{distance:5.0, sel:{model:lm}}},
             {stick:{radius:0.055, colorscheme:'grayCarbon', opacity:0.95}});

  const legend = [];
  spec.sidechains.forEach((sc) => {
    const m = v.addModel(sc.pdb, 'pdb');
    v.setStyle({model:m}, {stick:{radius:sc.radius, color:sc.color}});
    legend.push('<span style="color:' + sc.color + '">\\u25A0</span> ' + sc.label);
  });

  spec.measures.forEach((mz) => {
    v.addCylinder({start:{x:mz.a[0],y:mz.a[1],z:mz.a[2]},
                   end:{x:mz.b[0],y:mz.b[1],z:mz.b[2]},
                   radius:0.045, dashed:true, fromCap:1, toCap:1, color:mz.color});
    const mid = [0,1,2].map(i => (mz.a[i]+mz.b[i])/2);
    v.addLabel(mz.text, {position:{x:mid[0],y:mid[1],z:mid[2]},
                         backgroundColor:mz.color, backgroundOpacity:0.88,
                         fontColor:'#0d1117', fontSize:13, borderThickness:0,
                         inFront:true});
  });

  document.getElementById(keyId).innerHTML = legend.join('<br>');
  v.zoomTo({model:lm});
  v.zoom(D.zoom);
  v.render();
  return v;
}
build('p0','k0', D.panes[0]);
build('p1','k1', D.panes[1]);
window.__ready = true;
</script></body></html>
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="EPAS1")
    ap.add_argument("--chembl", default="CHEMBL4585668", help="belzutifan")
    ap.add_argument("--name", default="belzutifan")
    ap.add_argument("--resname", default="MET")
    ap.add_argument("--resid", type=int, default=252)
    ap.add_argument("--width", type=int, default=760)
    ap.add_argument("--height", type=int, default=620)
    ap.add_argument("--zoom", type=float, default=0.62)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    exp = ROOT / "exp" / args.gene
    rec_pdb = exp / "treatment" / "receptor.pdb"
    pose = exp / "flex" / "poses" / f"{args.chembl}.pdbqt"
    for p in (rec_pdb, pose, THREEDMOL):
        if not p.exists():
            print(f"[!] 없음: {p}")
            return 1

    lig, flex = parse_flex_pose(pose)
    key = f"{args.resname}{args.resid}"
    if key not in flex:
        print(f"[!] 포즈에 유연 잔기 {key} 가 없다. 있는 것: {list(flex)}")
        return 1

    moved = flex[key]
    orig = residue_atoms(rec_pdb, args.resname, args.resid)
    od = {n: p for n, p, _ in orig}
    common = [(n, p) for n, p, _ in moved if n in od]
    devs = [math.dist(od[n], p) for n, p in common]
    rms = math.sqrt(sum(d * d for d in devs) / len(devs))
    worst_d, worst_n = max(zip(devs, [n for n, _ in common]))

    # 접촉 거리 — 같은 리간드 포즈에 대해 두 측쇄 구조를 각각 잰다
    d_before, a_b, b_b = closest_pair(orig, lig)
    d_after, a_a, b_a = closest_pair(moved, lig)

    orig_sc = to_pdb(orig, args.resname, args.resid)
    moved_sc = to_pdb(moved, args.resname, args.resid)
    lig_pdb = to_pdb(lig, "LIG", 900, het=True)

    RED, GRN = "#f85149", "#3fb950"
    data = {
        "receptor": rec_pdb.read_text(),
        "ligand": lig_pdb,
        "zoom": args.zoom,
        "panes": [
            {
                "sidechains": [
                    {"pdb": orig_sc, "color": RED, "radius": 0.21,
                     "label": f"{key} 결정구조 (6CZW)"},
                ],
                "measures": [{"a": list(a_b), "b": list(b_b), "color": RED,
                              "text": f"{d_before:.2f} A"}],
            },
            {
                "sidechains": [
                    {"pdb": orig_sc, "color": RED, "radius": 0.10,
                     "label": f"{key} 결정구조"},
                    {"pdb": moved_sc, "color": GRN, "radius": 0.21,
                     "label": f"{key} 유연 도킹 후"},
                ],
                "measures": [{"a": list(a_a), "b": list(b_a), "color": GRN,
                              "text": f"{d_after:.2f} A"}],
            },
        ],
    }

    html = (
        PAGE.replace("__LIB__", THREEDMOL.read_text())
        .replace("__DATA__", json.dumps(data))
        .replace("__TOTW__", str(args.width * 2 + 1))
        .replace("__W__", str(args.width))
        .replace("__H__", str(args.height))
        .replace("__T0__", f"리지드 — {key} 고정")
        .replace("__S0__", f"{args.name} 가 들어오지 못하는 자리")
        .replace("__C0__", f"{key} ↔ {args.name}  {d_before:.2f} A\n"
                           f"탄소-탄소 반데르발스 접촉 3.4~4.0 A — 충돌")
        .replace("__T1__", f"유연 측쇄 — {key} 자유")
        .replace("__S1__", "측쇄가 회전해 충돌이 해소됨")
        .replace("__C1__", f"{key} ↔ {args.name}  {d_after:.2f} A  "
                           f"({d_after - d_before:+.2f})\n"
                           f"측쇄 RMSD {rms:.2f} A — {worst_n} {worst_d:.2f} A 이동\n"
                           f"PHE254 0.15 A, HIS293 0.37 A — 함께 풀었으나 안 움직임")
    )

    outdir = exp / "flex"
    tmp = outdir / "_m252.html"
    tmp.write_text(html)
    out = Path(args.out) if args.out else outdir / f"{key.lower()}_before_after.png"

    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        b = pw.chromium.launch(args=["--no-sandbox", "--use-gl=swiftshader",
                                     "--enable-unsafe-swiftshader"])
        pg = b.new_page(viewport={"width": args.width * 2 + 1, "height": args.height},
                        device_scale_factor=2)
        pg.goto(tmp.as_uri())
        pg.wait_for_function("window.__ready === true", timeout=60000)
        pg.wait_for_timeout(2500)
        pg.screenshot(path=str(out))
        b.close()
    tmp.unlink(missing_ok=True)

    print(f"[*] {key} 측쇄 이동")
    print(f"    RMSD {rms:.2f} A, 최대 {worst_n} {worst_d:.2f} A")
    print(f"    {args.name} 접촉 {d_before:.2f} -> {d_after:.2f} A ({d_after - d_before:+.2f})")
    print(f"[+] {out}  ({out.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
