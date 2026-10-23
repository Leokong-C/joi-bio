# Joi_bio — EPAS1 (HIF-2 alpha) PAS-B
#
# PyMOL 에서:  pymol joi_bio.pml
# 또는 PyMOL 안에서:  @joi_bio.pml
#
# 장면 세 개가 들어 있다. F1~F3 또는 scene 메뉴로 전환.

reinitialize
bg_color white
set ray_opaque_background, 0

# ---------- 구조 ----------
load receptor_rigid.pdb,     rigid
load receptor_m252open.pdb,  m252open
load poses_rigid.sdf,        poses_rigid

# 유연 측쇄 도킹 결과 (벤치마크 3개)
load flex_belzutifan_lig.pdb,  flex_belz_lig
load flex_belzutifan_sc.pdb,   flex_belz_sc
load flex_PT2399_lig.pdb,      flex_pt2399_lig
load flex_PT2385_lig.pdb,      flex_pt2385_lig

# ---------- 기본 표현 ----------
hide everything
show cartoon, rigid or m252open
color grey80, rigid
color palecyan, m252open
set cartoon_transparency, 0.45, rigid
set cartoon_transparency, 0.45, m252open

# 문헌이 지목한 세 잔기 (Ren et al. 2022, Mol Pharmacol 102:240)
select key_res_rigid, rigid and resi 252+254+293
select key_res_open,  m252open and resi 252+254+293
show sticks, key_res_rigid or key_res_open
color grey50, key_res_rigid and elem C
color deepteal, key_res_open and elem C
util.cnc("key_res_rigid or key_res_open")

# ---------- 장면 1: 리지드 도킹 포즈 ----------
disable m252open
disable flex_belz_lig
disable flex_belz_sc
disable flex_pt2399_lig
disable flex_pt2385_lig
show sticks, poses_rigid
color yellow, poses_rigid and elem C
util.cnc("poses_rigid")
orient rigid and resi 252+254+293
zoom rigid and resi 252+254+293, 8
scene overview, store, message="리지드 도킹 — 벤치마크 3개 + 상위 히트 20개 (노란색). 회색 막대가 M252/F254/H293"

# ---------- 장면 2: M252 측쇄 뒤집힘 ----------
disable poses_rigid
enable m252open
show sticks, flex_belz_lig
color orange, flex_belz_lig and elem C
util.cnc("flex_belz_lig")
enable flex_belz_lig
distance d_rigid, rigid and resi 252 and name CE, flex_belz_lig
color red, d_rigid
hide labels, d_rigid
orient rigid and resi 252
zoom rigid and resi 252, 6
scene m252_flip, store, message="M252 뒤집힘 — 회색이 원본 6CZW, 청록이 유연 도킹 후. 측쇄 RMSD 2.37 A, CE 원자 3.96 A 이동. 접촉 거리 1.53 -> 3.21 A"

# ---------- 장면 3: 유연 측쇄 도킹 결과 ----------
disable rigid
enable m252open
show sticks, flex_belz_sc
color salmon, flex_belz_sc and elem C
util.cnc("flex_belz_sc")
enable flex_belz_sc
enable flex_pt2399_lig
enable flex_pt2385_lig
show sticks, flex_pt2399_lig or flex_pt2385_lig
color marine, flex_pt2399_lig and elem C
color purple, flex_pt2385_lig and elem C
util.cnc("flex_pt2399_lig or flex_pt2385_lig")
orient m252open and resi 252+254+293
zoom m252open and resi 252+254+293, 8
scene flex_result, store, message="유연 측쇄 도킹 — belzutifan(주황) PT2399(파랑) PT2385(보라). 연어색이 움직인 측쇄. 세 잔기를 풀었으나 실제로 움직인 것은 M252 하나"

# ---------- 마무리 ----------
scene overview, recall
set stick_radius, 0.14
set cartoon_side_chain_helper, on
print ""
print "=== Joi_bio EPAS1 ==="
print "장면: overview / m252_flip / flex_result   (F1~F3 또는 scene 메뉴)"
print ""
print "유연 측쇄 1단계 결과 (실행 38973161302, 표본 53개)"
print "  AUC      0.7667 -> 0.8478"
print "  BEDROC   0.9656 -> 0.9978"
print "  벤치마크 순위 중앙값 +9 계단"
print "  판정     INDUCED_FIT_SUPPORTED"
print ""
print "움직인 잔기는 M252 하나뿐이다 (측쇄 RMSD 2.37 A)."
print "PHE254 는 0.16 A, HIS293 은 0.40 A 로 사실상 제자리."
print ""
