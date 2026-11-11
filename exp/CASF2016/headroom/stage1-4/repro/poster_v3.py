
import matplotlib as mpl, matplotlib.pyplot as plt, numpy as np, pandas as pd, textwrap, os, json
mpl.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.edgecolor":"#cccccc",
    "axes.linewidth":0.8,"grid.color":"#eeeeee","xtick.color":"#333333","ytick.color":"#333333",
    "text.color":"#1a1a1a","font.family":"Malgun Gothic","axes.unicode_minus":False,"font.size":9})
N = json.load(open("poster_numbers.json"))
mas = pd.read_csv("casf2016_predictor_matrix.csv")
imp = pd.read_csv("casf2016_null_feature_importance.csv", index_col=0)
dec = pd.read_csv("casf2016_null_decomposition_corrected.csv")
seeds = pd.read_csv("casf2016_seed_repeats.csv")
srt = seeds.groupby("seed").apply(lambda g: 100*(g.rmsd_top1<=2).mean(), include_groups=False)
_r = dec[(dec.predictor=="Vina") & (dec.controlled_for=="ligand size")].iloc[0]
R2_size, R2_vina, U_vina = float(_r.R2_control_alone), float(_r.R2_predictor_alone), float(_r.unique_R2_added)
R2_both = R2_size + U_vina; U_size = R2_both - R2_vina; SHARED = R2_both - U_vina - U_size

V="#1f6feb"; VD="#8e44ad"; NU="#c48f1b"; R="#d1495b"; G="#8a8f98"; ML="#0b6b5e"; LM="#8a4b9e"; INK="#1a1a1a"
WRAP_T, WRAP_F = 54, 112

STAGES = [
 dict(tag="①", col=NU, title="널 모델 기준선 — 아무 물리도 모르는 예측자",
   goal="결합을 전혀 모르는 예측자가 어디까지 가는지 먼저 정한다.",
   res=(f"리간드 중원자 수 하나만으로 log$K_a$와 $R$ = +{N['R_size']:.3f}, 타깃 내 순위 "
        f"상관 ρ = +{N['rho_size']:.3f}. 이 ρ는 이 포스터의 어떤 스코어링 함수보다도 높다."),
   caveat="결정 해상도도 log$K_a$와 약한 상관(R = +0.173)을 보인다 — 교란 요인."),
 dict(tag="②", col=V, title="Vina 재도킹 — 세 가지 능력을 분리한다",
   goal="포즈 복원 · 세기 예측 · 같은 타깃 내 순위를 따로 측정한다.",
   res=(f"top-1이 2 Å 이내 {N['top1']:.1f} % (중앙값 {N['med']:.2f} Å), 반환 모드 중 최선은 "
        f"{N['best']:.1f} %. 물리 타당성 {N['pbv']:.1f} %, 타당성 AND 2 Å는 {N['pbb']:.1f} %$^{{7,17}}$. "
        f"scoring $R$ = +{N['R_vina']:.3f}, ranking ρ = +{N['rho_vina']:.3f}."),
   caveat=f"ranking ρ가 크기 널 모델(+{N['rho_size']:.3f})보다 낮다 — 순위에서는 널에 진다."),
 dict(tag="③", col=NU, title="크기 효과를 걷어내면 무엇이 남는가",
   goal="'Vina가 널을 겨우 넘는다'가 'Vina는 크기만 본다'를 뜻하는지 본다.",
   res=(f"Vina 점수와 리간드 크기의 상관이 $R$ = +0.487. 크기를 통제한 편상관은 "
        f"+{N['pr_size_vina']:.3f}" r" ($p$ = 3$\times$10$^{-9}$)로 여전히 유의하다. "
        f"분산분해: 크기 고유 {U_size:.3f} · 공유 {SHARED:.3f} · Vina 고유 {U_vina:.3f}."),
   caveat="타깃 내부에서는 평균 R 0.501 → 편상관 0.224로 떨어진다."),
 dict(tag="④", col=VD, title="스코어링 함수를 바꾸면 해결되는가",
   goal="같은 포즈를 Vinardo로 재채점해 함수 선택의 문제인지 본다.$^{6}$",
   res=(f"Vinardo $R$ = +{N['R_vinardo']:.3f} — Vina(+{N['R_vina']:.3f})보다도 크기 널"
        f"(+{N['R_size']:.3f})보다도 낮다. 포즈 순위는 {N['reord']:.1f} %에서 바뀌지만 "
        f"docking power는 {N['top1']:.1f} % → {N['top1_vd']:.1f} %뿐."),
   caveat=f"두 함수의 점수 상관이 $R$ = 0.942 — 사실상 같은 신호를 본다."),
 dict(tag="⑤", col=V, title="탐색이 부족한가, 순위가 틀리는가",
   goal="탐색량을 4배로 늘려 top-1과 '최선 모드'의 간격이 좁혀지는지 본다.",
   res=(f"exhaustiveness 8 → 32 ({N['cpu']:.0f} CPU-시간): top-1 성공률 {srt.mean():.1f} % → "
        f"{N['ex32']:.1f} %. 시드 3회 반복의 폭은 {N['s_lo']:.1f}–{N['s_hi']:.1f} % (SD {N['s_sd']:.2f} pt) — "
        f"4배 탐색의 이득이 잡음에 묻힌다."),
   caveat=f"복합체 {N['flip']}건(9.5 %)은 시드만 바꿔도 2 Å 판정이 뒤집힌다.$^{{16}}$"),
 dict(tag="⑥", col=NU, title="다변량 널 모델 — 사소한 기술자를 전부 동원하면",
   goal="결합을 모르는 기술자 14개를 전부 쓰면 어디까지 가는지 본다.",
   res=(f"표적으로 묶은 GroupKFold(10)에서 랜덤포레스트 $R$ = +{N['R_rf']:.3f}, 릿지 +{N['R_ridge']:.3f} — "
        f"Vina와 Vinardo를 모두 넘는다. 최대 기여 특징은 리간드 6 Å 안의 단백질 원자 개수"
        f"({N['imp_top']:.1f} %) — 화학이 아니라 주머니 밀도다."),
   caveat="널이 이긴다는 것은 '널이 좋다'가 아니라 '그 지표가 증거가 못 된다'는 뜻이다."),
 dict(tag="⑦", col=ML, title="기계학습 스코어링 함수는 다른가",
   goal="PDBbind refined에서 코어셋을 뺀 데이터로 학습해 적용한다.$^{10,11,12}$",
   res=(f"RF-Score v2 $R$ = +{N['R_rfs2']:.3f} (v1 +{N['R_rfs1']:.3f}), 잔차 SD {N['sd_rfs2']:.2f} vs Vina "
        f"{N['sd_vina']:.2f}. 널을 제거하고 남는 고유 기여가 Δ$R^2$ = {N['u_rfs2']:.3f}로, Vina "
        f"{N['u_vina']:.3f}·Vinardo {N['u_vd']:.3f}와 차원이 다르다."),
   caveat=f"타깃 내 순위는 ρ = +{N['rho_rfs2']:.3f}로 중원자 수 세기(+{N['rho_size']:.3f})에 여전히 진다."),
 dict(tag="⑧", col=LM, title="대규모 언어모델 — 제로샷과 기권",
   goal="SMILES와 단백질 이름만 주고, 근거가 없으면 기권하게 한다.$^{17}$",
   res=(f"리간드만 주면 $R$ = +{N['llm_A']:.3f}이고 {N['ab_A']:.0f} % 기권한다. 단백질 이름을 주면 "
        f"+{N['llm_B']:.3f}(기권 {N['ab_B']:.0f} %), PDB ID까지 주면 +{N['llm_C']:.3f}(기권 {N['ab_C']:.0f} %) — "
        f"Vina와 널을 앞서고 RF-Score v2에는 못 미친다."),
   caveat=(f"오염 검증 음성 — '안다'고 보고한 쪽이 오히려 더 틀린다($R$ {N['rec_yes']:.3f} vs {N['rec_no']:.3f}). "
           "PDB ID가 사는 것은 정확도가 아니라 자신감이다.")),
]

foot_raw = [
 "1  Su M. et al. Comparative Assessment of Scoring Functions: The CASF-2016 Update. J Chem Inf Model 2018. doi:10.1021/acs.jcim.8b00545",
 "2  Su M. et al. Assessing protein-ligand interaction scoring functions with the CASF-2013 benchmark. Nat Protoc 2018. doi:10.1038/nprot.2017.114",
 "3  Wang R. et al. The PDBbind Database. J Med Chem 2004. doi:10.1021/jm030580l",
 "4  Trott O., Olson A. J. AutoDock Vina. J Comput Chem 2010. doi:10.1002/jcc.21334",
 "5  Eberhardt J. et al. AutoDock Vina 1.2.0. J Chem Inf Model 2021. doi:10.1021/acs.jcim.1c00203",
 "6  Quiroga R., Villarreal M. A. Vinardo. PLOS ONE 2016. doi:10.1371/journal.pone.0155183",
 "7  Buttenschoen M. et al. PoseBusters. Chem Sci 2024. doi:10.1039/d3sc04185a",
 "8  Meli R., Biggin P. C. spyrmsd: symmetry-corrected RMSD calculations in Python. J Cheminform 2020. doi:10.1186/s13321-020-00455-2",
 "9  O'Boyle N. M. et al. Open Babel: An open chemical toolbox. J Cheminform 2011. doi:10.1186/1758-2946-3-33",
 "10 Ballester P. J., Mitchell J. B. O. A machine learning approach to predicting protein-ligand binding affinity. Bioinformatics 2010. doi:10.1093/bioinformatics/btq112",
 "11 Ballester P. J., Schreyer A., Blundell T. L. Does a More Precise Chemical Description Lead to More Accurate Prediction? J Chem Inf Model 2014. doi:10.1021/ci500091r",
 "12 Wojcikowski M., Zielenkiewicz P., Siedlecki P. Open Drug Discovery Toolkit (ODDT). J Cheminform 2015. doi:10.1186/s13321-015-0078-2",
 "13 방법: 수용체는 PDB의 ATOM 레코드만 남겨 obabel -xr -p 7.4로 PDBQT 변환. 리간드는 RCSB ModelServer가 주는 결합차수 포함 SDF를 Meeko로 PDBQT 변환.",
 "14 n = 285: 1a30(펩타이드 리간드)을 포함하고 3syr의 리간드를 보조인자 LLP에서 실제 리간드 GPK로 정정했다. 리간드 식별은 RCSB 결합친화도 주석이 있는 50개에서 49/50 일치로 검증.",
 "15 정정: mmCIF HETATM 중원자 집계는 대체위치(altloc)를 중복으로 센다(3gv9·2xnb·2yfe·2al5는 정확히 2배, 3prs는 50 대 53). SDF 기준으로 세면 크기 널 모델 R이 +0.469에서 +0.482로 바뀐다.",
 "16 재현성 점검: 같은 설정 재실행에서 284개 중 250개가 소수점 3자리까지 일치, 성공 판정 뒤집힘 0건. 차이의 원인은 Open Babel의 프로세스 간 비결정적 극성 수소 배치로 확정.",
 "17 RF-Score는 ODDT 0.4로 PDBbind v2016 refined set에서 코어셋을 제외해 학습(train_set='refined' & ~core, 코드로 확인). 언어모델은 claude-sonnet-5, 사고과정 비활성, 도구 호출 강제.",
 "18 PoseBusters 수치는 n = 284 — 정정한 2개 복합체는 재검사하지 않았다. 상자 크기 민감도는 층화 추출한 60개 부분집합에서 측정했다.",
 "모든 DOI는 Crossref에서 해석 확인. 문헌 메타데이터는 PubMed 조회.",
]

def wrap(s, w): return textwrap.wrap(s, w, break_long_words=False)
LH_G, LH_R, LH_C = 0.175, 0.183, 0.175
for st in STAGES:
    st["g"] = wrap(st["goal"], WRAP_T); st["r"] = wrap(st["res"], WRAP_T); st["c"] = wrap(st["caveat"], WRAP_T)
    _bodyh = 0.50 + LH_G*len(st["g"]) + 0.10 + LH_R*len(st["r"]) + 0.10 + LH_C*len(st["c"]) + 0.24
    st["bh"] = max(_bodyh, 1.62)

SUMMARY = (
 f"① 포즈 생성은 병목이 아니다 — 올바른 포즈는 {N['best']:.1f} %의 복합체에서 만들어지지만 점수가 1위로 올리는 것은 {N['top1']:.1f} %뿐이고, 탐색량을 4배({N['cpu']:.0f} CPU-시간)로 써도 간격은 줄지 않는다.\n"
 f"② 물리 기반 스코어링은 사소한 기술자에 거의 흡수된다 — 널 모델을 제거하고 남는 고유 기여가 Vina Δ$R^2$ = {N['u_vina']:.3f}, Vinardo {N['u_vd']:.3f}. 함수를 Vinardo로 바꾸는 것은 해결책이 아니다.\n"
 f"③ 기계학습과 언어모델은 다르다 — RF-Score v2는 같은 널을 제거하고도 Δ$R^2$ = {N['u_rfs2']:.3f}를 더하고, 제로샷 언어모델도 Vina를 앞선다(+{N['llm_C']:.3f}).\n"
 f"④ 그런데 타깃 내 순위에서는 여전히 '중원자 수 세기'(ρ = +{N['rho_size']:.3f})가 1위다.  →  남은 문제는 세기 예측이 아니라 같은 표적 안에서의 순위다.")
CAVEATS = (
 f"· 단일 실행 수치를 소수점까지 비교하면 안 된다 — 시드만 바꿔도 복합체 {N['flip']}건(9.5 %)의 2 Å 판정이 뒤집히고 성공률 표준편차는 {N['s_sd']:.2f} pt다.$^{{16}}$\n"
 f"· 탐색 상자 크기는 시드 잡음보다 훨씬 크게 작용한다 — 여유 +10 Å → +15 Å만으로 top-1 성공률이 {N['bx10']:.1f} % → {N['bx15']:.1f} %, 최선 모드는 {N['bx_best5']:.1f} % → {N['bx_best15']:.1f} %로 떨어진다.$^{{18}}$\n"
 f"· 비트 단위 재현에는 준비된 수용체 PDBQT 보관이 필요하다 — Open Babel의 극성 수소 배치가 프로세스 간 비결정적이다.$^{{9,16}}$\n"
 f"· 리간드 식별이 틀리면 거짓 실패가 생긴다 — 3syr을 보조인자 LLP로 잘못 잡으면 top-1 RMSD 4.41 Å(실패)이지만 올바른 리간드 GPK로는 0.27 Å(성공)이다.$^{{14}}$\n"
 f"· 널 모델 비교는 CASF 코어셋의 구성(타깃당 리간드 5개, 넓은 크기 분포)에 의존한다 — 다른 세트에서는 크기 효과가 달라질 수 있다.\n"
 f"· 공동접힘(co-folding) 모델(DiffDock·Boltz·AlphaFold3)과 자유에너지 섭동, 스크리닝 능력(screening power)은 이 포스터의 범위 밖이다.")
foot = []
for f_ in foot_raw: foot += wrap(f_, 148) if len(f_) > 148 else [f_]

HEAD = 0.42+0.50+0.33+0.28+0.42
BAND = 1.00+0.18
STACK = sum(s["bh"] for s in STAGES) + 0.105*(len(STAGES)-1) + 0.12
SUMH = 0.34 + 0.225*len(SUMMARY.split("\n")) + 0.18
CAVH = 0.21 + 0.165*len(CAVEATS.split("\n")) + 0.20
FOOTH = 0.20 + 0.152*len(foot)
W = 11.7
H = HEAD + BAND + STACK + SUMH + 0.16 + CAVH + FOOTH + 0.34
print("computed height: %.2f in  (head %.2f band %.2f stack %.2f summary %.2f caveats %.2f foot %.2f)"
      % (H, HEAD, BAND, STACK, SUMH, CAVH, FOOTH))

_probe_fig = plt.figure(figsize=(W, H)); _probe_r = _probe_fig.canvas.get_renderer()
def _wpx(s, fontsize, weight="normal"):
    t = _probe_fig.text(0, 0, s, fontsize=fontsize, weight=weight)
    w = t.get_window_extent(_probe_r).width; t.remove(); return w
def wrap_fit(s, x0, x1, fontsize, weight="normal"):
    """Greedy wrap so every line fits between figure fractions x0 and x1."""
    lim = (x1 - x0) * _probe_fig.bbox.width
    out, cur = [], ""
    for word in s.split(" "):
        cand = word if not cur else cur + " " + word
        if _wpx(cand, fontsize, weight) <= lim or not cur:
            cur = cand
        else:
            out.append(cur); cur = word
    if cur: out.append(cur)
    return out
def wrap_block(text, x0, x1, fontsize, weight="normal"):
    lines = []
    for para in text.split("\n"): lines += wrap_fit(para, x0, x1, fontsize, weight)
    return lines

for st in STAGES:
    st["g"]=wrap_fit(st["goal"], 0.133, 0.600, 9.2)
    st["r"]=wrap_fit(st["res"],  0.133, 0.600, 9.6)
    st["c"]=wrap_fit(st["caveat"],0.148, 0.600, 9.0)
    st["bh"]=max(0.50+LH_G*len(st["g"])+0.10+LH_R*len(st["r"])+0.10+LH_C*len(st["c"])+0.24, 1.62)
SUM_L  = wrap_block(SUMMARY, 0.055, 0.945, 9.8)
CAV_L  = wrap_block(CAVEATS, 0.040, 0.958, 8.9)
FOOT_L = []
for f_ in foot_raw: FOOT_L += wrap_fit(f_, 0.040, 0.958, 7.3)
plt.close(_probe_fig)
STACK = sum(s["bh"] for s in STAGES) + 0.105*(len(STAGES)-1) + 0.12
SUMH  = 0.38 + 0.225*len(SUM_L) + 0.20
CAVH  = 0.21 + 0.165*len(CAV_L) + 0.20
FOOTH = 0.20 + 0.152*len(FOOT_L)
H = HEAD + BAND + STACK + SUMH + 0.16 + CAVH + FOOTH + 0.34
print("final height %.2f in | summary %d lines, caveats %d, footnotes %d" % (H, len(SUM_L), len(CAV_L), len(FOOT_L)))

def fy(t): return 1.0 - t/H
def fh(i): return i/H
plt.close("all")
fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor("white")
def T(x, t, s, **kw):
    kw.setdefault("ha","left"); kw.setdefault("va","top"); kw.setdefault("color",INK); kw.setdefault("zorder",5)
    return fig.text(x, fy(t), s, **kw)
def rect(x, t_bottom, w, h_in, **kw):
    fig.patches.append(mpl.patches.Rectangle((x, fy(t_bottom)), w, fh(h_in), transform=fig.transFigure, **kw))

y = 0.42
T(0.04, y, "도킹 프로그램은 자리를 맞히는가, 세기를 맞히는가", fontsize=25, weight="bold"); y += 0.50
fig.add_artist(mpl.lines.Line2D([0.04,0.96],[fy(y-0.09)]*2, color=INK, lw=2.4, zorder=5))
T(0.04, y, f"CASF-2016 코어셋 {N['n']}개 복합체에서 물리 기반·기계학습·언어모델 스코어링을 같은 자로 재다$^{{1,4,5}}$",
  fontsize=13.5, color="#444"); y += 0.33
T(0.04, y, "핵심 질문 — 스코어링 함수가 친화도와 상관이 있다는 사실이, 그 함수가 결합을 모형화한다는 증거가 되는가?",
  fontsize=12, color=V); y += 0.28
T(0.04, y, "답 — 물리 기반 함수에게는 거의 증거가 되지 못한다. 결합을 전혀 모르는 기술자로 만든 널 모델이 Vina와 Vinardo를 둘 다 이긴다.",
  fontsize=12, weight="bold", color=R); y += 0.42

rect(0.04, y+1.00, 0.92, 1.00, facecolor="#f2f4f7", edgecolor="#d0d5dd", lw=.9, zorder=0)
T(0.052, y+0.07, "설계", fontsize=11.5, weight="bold")
T(0.105, y+0.07,
 f"CASF-2016 코어셋$^{{1,3}}$ · 57개 타깃 × 리간드 5개 = {N['n']}$^{{14}}$ · 실험 log$K_a$ {N['lo']}–{N['hi']}\n"
 f"수용체는 ATOM 레코드만 → Open Babel pH 7.4 PDBQT$^{{9}}$ · 리간드는 RCSB ModelServer SDF → Meeko PDBQT · 상자 = 리간드 범위 + 10 Å$^{{13}}$\n"
 f"재도킹 exhaustiveness 8, num_modes 9, seed 42, 작업당 1 CPU · RMSD는 sPyRMSD 대칭 보정$^{{8}}$ · 지표 정의는 CASF 프로토콜$^{{2}}$\n"
 f"비교 대상 7종 — 물리 기반 2(Vina, Vinardo) · 사소한 널 3(중원자 수, 릿지, 랜덤포레스트) · 기계학습 2(RF-Score v1·v2) · 언어모델 1(제로샷·기권)",
 fontsize=9.3, color="#333", linespacing=1.52)
y += 1.18

panel_axes = []
for st in STAGES:
    top = y; bh = st["bh"]; bot = top + bh
    rect(0.04, bot, 0.92, bh, facecolor="white", edgecolor="#d0d5dd", lw=.9, zorder=0)
    rect(0.04, bot, 0.0065, bh, facecolor=st["col"], edgecolor="none", zorder=1)
    T(0.058, top+0.11, st["tag"], fontsize=16.5, weight="bold", color=st["col"])
    T(0.090, top+0.13, st["title"], fontsize=13.5, weight="bold")
    yy = top + 0.50
    T(0.090, yy, "목적", fontsize=9, weight="bold", color="#666")
    T(0.133, yy, "\n".join(st["g"]), fontsize=9.2, color="#444", linespacing=1.5)
    yy += LH_G*len(st["g"]) + 0.10
    T(0.090, yy, "결과", fontsize=9, weight="bold", color="#666")
    T(0.133, yy, "\n".join(st["r"]), fontsize=9.6, linespacing=1.5)
    yy += LH_R*len(st["r"]) + 0.10
    T(0.090, yy, "유의점", fontsize=9, weight="bold", color=R)
    T(0.148, yy, "\n".join(st["c"]), fontsize=9.0, color=R, linespacing=1.5)
    ax = fig.add_axes([0.655, fy(bot-0.26), 0.285, fh(bh-0.58)])
    ax.set_zorder(4); ax.patch.set_alpha(0); ax.tick_params(labelsize=7.5)
    [ax.spines[s].set_visible(False) for s in ("top","right")]
    panel_axes.append(ax)
    y = bot + 0.105
y += 0.12

# ── panels ───────────────────────────────────────────────────────────────────
ax = panel_axes[0]                                   # ① null scatter
ax.plot(mas.n_heavy, mas.logKa, "o", ms=2.6, color=G, alpha=.7, mew=0)
b_, a_ = np.polyfit(mas.n_heavy, mas.logKa, 1)
xs = np.linspace(mas.n_heavy.min(), mas.n_heavy.max(), 40)
ax.plot(xs, a_+b_*xs, color=NU, lw=1.6)
ax.set_xlabel("리간드 중원자 수", fontsize=8); ax.set_ylabel("log$K_a$", fontsize=8)
ax.text(.03,.97, f"널 모델  $R$ = +{N['R_size']:.3f}", transform=ax.transAxes, va="top", fontsize=8.5, color=NU)

ax = panel_axes[1]                                   # ② docking power
vals=[N["top1"], N["best"], N["pbb"]]
ax.bar(range(3), vals, color=[V,G,"#2e8b57"], width=.6, edgecolor="white", lw=.7)
for i_,v in enumerate(vals): ax.text(i_, v+2.5, f"{v:.1f}%", ha="center", fontsize=8.5)
ax.set_xticks(range(3)); ax.set_xticklabels(["top-1\n2 Å","최선 모드\n2 Å","타당성\nAND 2 Å"], fontsize=7.5)
ax.set_ylim(0,108); ax.set_ylabel("복합체 비율 (%)", fontsize=8)

ax = panel_axes[2]                                   # ③ variance decomposition
left=0
for lab,val,c in [("크기 고유",U_size,NU),("공유",SHARED,G),("Vina 고유",U_vina,V)]:
    ax.barh([0],[val], left=left, color=c, height=.40, edgecolor="white", lw=.8)
    ax.text(left+val/2, 0, f"{val:.3f}", ha="center", va="center", fontsize=8, color="white")
    ax.text(left+val/2, 0.30, lab, ha="center", fontsize=8, color=c); left += val
ax.set_xlim(0,0.37); ax.set_ylim(-0.42,0.62); ax.set_yticks([])
ax.set_xlabel("log$K_a$ 설명 분산 ($R^2$)", fontsize=8); ax.spines["left"].set_visible(False)

ax = panel_axes[3]                                   # ④ Vina vs Vinardo vs size
pe=[N["R_vina"],N["R_vinardo"],N["R_size"]]; pa=[N["pr_size_vina"],0.300,np.nan]
xs3=np.arange(3); cc=[V,VD,NU]
ax.bar(xs3-0.17, pe, width=.32, color=cc, edgecolor="white", lw=.7)
ax.bar(xs3+0.17, [0 if p!=p else p for p in pa], width=.32, color=cc, alpha=.40, hatch="///", edgecolor="white", lw=.7)
for i_,(p1,p2) in enumerate(zip(pe,pa)):
    ax.text(i_-0.17, p1+0.014, f"{p1:.3f}", ha="center", fontsize=7.5)
    if p2==p2: ax.text(i_+0.17, p2+0.014, f"{p2:.3f}", ha="center", fontsize=7.5)
ax.set_xticks(xs3); ax.set_xticklabels(["Vina","Vinardo","리간드 크기"], fontsize=8)
ax.set_ylabel("log$K_a$와의 Pearson $R$", fontsize=8); ax.set_ylim(0,0.62)
ax.text(.03,.97,"빗금 = 크기 제거 후", transform=ax.transAxes, va="top", fontsize=7.5)

ax = panel_axes[4]                                   # ⑤ exhaustiveness
ax.axhspan(N["s_lo"], N["s_hi"], color=V, alpha=.15, zorder=0)
ax.plot([8]*3, srt.values, "o", ms=5, color=V, mew=0)
ax.plot([32],[N["ex32"]], "o", ms=5, color=V, mew=0)
ax.plot([8,32],[srt.mean(), N["ex32"]], "-", color=V, lw=1.3)
ax.set_xscale("log", base=2); ax.set_xticks([8,32]); ax.set_xticklabels(["8","32"], fontsize=8)
ax.set_xlim(6,44); ax.set_ylim(63.5,68.5)
ax.set_xlabel("exhaustiveness (탐색량)", fontsize=8); ax.set_ylabel("top-1 2 Å 성공률 (%)", fontsize=8)
ax.text(.03,.96, f"음영 = 시드 3회 폭\n4배 탐색 → +{N['ex32']-srt.mean():.2f} pt", transform=ax.transAxes,
        va="top", fontsize=8, color=V)

ax = panel_axes[5]                                   # ⑥ feature importance
top6 = imp.iloc[:6,0][::-1]
nice = {"c60":"6 Å 내 단백질 원자","logp":"cLogP","frac_csp3":"sp$^3$ 탄소 비율","mw":"분자량",
        "rgyr":"회전반경","tpsa":"극성 표면적","c45_per_heavy":"중원자당 접촉","c45":"4.5 Å 내 원자"}
ax.barh(np.arange(6), top6.values*100, color=[NU if i<5 else "#9a7318" for i in range(6)], height=.66)
ax.set_yticks(np.arange(6)); ax.set_yticklabels([nice.get(i,i) for i in top6.index], fontsize=7.5)
ax.set_xlabel("랜덤포레스트 특징 중요도 (%)", fontsize=8); ax.set_xlim(0,36)
ax.text(top6.values[-1]*100-1, 5, f"{top6.values[-1]*100:.1f}", ha="right", va="center",
        fontsize=8, color="white", weight="bold")

ax = panel_axes[6]                                   # ⑦ RF-Score predicted vs experimental
ax.plot([2,12],[2,12], ls="--", lw=1.0, color=G, zorder=1)
ax.plot(mas.logKa, mas.rfscore_v2, "o", ms=2.6, color=ML, alpha=.70, mew=0, zorder=3)
ax.set_xlim(1.7,12.3); ax.set_ylim(1.7,12.3)
ax.set_xlabel("실험 log$K_a$", fontsize=8); ax.set_ylabel("RF-Score v2 예측 p$K_d$", fontsize=8)
ax.text(.03,.97, f"$R$ = +{N['R_rfs2']:.3f}\n잔차 SD {N['sd_rfs2']:.2f}", transform=ax.transAxes,
        va="top", fontsize=8.5, color=ML)

ax = panel_axes[7]                                   # ⑧ LLM conditions
rr=[N["llm_A"],N["llm_B"],N["llm_C"]]; ab=[N["ab_A"],N["ab_B"],N["ab_C"]]
ax.bar(range(3), rr, color=["#c9a3d6",LM,LM], width=.58, edgecolor="white", lw=.7)
for i_,(v,a2) in enumerate(zip(rr,ab)):
    ax.text(i_, v+0.015, f"{v:.3f}", ha="center", fontsize=8.5, color=LM if v>0.5 else "#9a7aa8")
    ax.text(i_, 0.025, f"기권 {a2:.0f}%", ha="center", fontsize=7.5, color="white")
ax.axhline(N["R_vina"], color=V, lw=1.1, ls=":", zorder=4)
ax.text(2.46, N["R_vina"], "Vina", fontsize=7.5, color=V, va="center", ha="right")
ax.set_xticks(range(3)); ax.set_xticklabels(["리간드만","+단백질명","+PDB ID"], fontsize=8)
ax.set_ylabel("log$K_a$와의 Pearson $R$", fontsize=8); ax.set_ylim(0,0.80)


rect(0.04, y+SUMH-0.20, 0.92, SUMH-0.20, facecolor="#eef3fb", edgecolor=V, lw=1.1, zorder=0)
T(0.055, y+0.09, "종합", fontsize=13.5, weight="bold", color=V)
T(0.055, y+0.40, "\n".join(SUM_L), fontsize=9.8, linespacing=1.62)
y += SUMH + 0.16
T(0.04, y, "유의점 (전체)", fontsize=10.5, weight="bold", color=R); y += 0.21
T(0.04, y, "\n".join(CAV_L), fontsize=8.9, color="#333", linespacing=1.60)
y += CAVH - 0.21
T(0.04, y, "각주", fontsize=10.5, weight="bold"); y += 0.20
for j, f_ in enumerate(FOOT_L): T(0.04, y + 0.152*j, f_, fontsize=7.3, color="#444")
end = y + 0.152*len(FOOT_L)
fig.savefig("casf2016_poster.png", dpi=220, facecolor="white")
r_ = fig.canvas.get_renderer()
tx = [(t, t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
oob = [t.get_text()[:26] for t,b in tx if not (fig.bbox.contains(b.x0,b.y0) and fig.bbox.contains(b.x1,b.y1))]
print("saved %.2f MB | %.1f x %.1f in | ends %.2f in (margin %.2f) | out of page: %s"
      % (os.path.getsize("casf2016_poster.png")/1e6, W, H, end, H-end, oob[:6] or "none"))
