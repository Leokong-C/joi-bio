Joi_bio — PyMOL 묶음 (EPAS1 / HIF-2 alpha PAS-B)
생성 2026-10-23

여는 법
  pymol joi_bio.pml
  또는 PyMOL 안에서  @joi_bio.pml

파일
  receptor_rigid.pdb        6CZW chain A, 880 원자. 리지드 도킹에 쓴 원본 수용체
  receptor_m252open.pdb     같은 구조에서 M252 측쇄만 유연 도킹 결과로 바꾼 것
                            880 원자 보존, MET252 의 CG/SD/CE 만 좌표가 다르다
  poses_rigid.sdf           리지드 도킹 포즈 23개 (벤치마크 3 + 상위 히트 20)
                            각 분자에 vina_score, group, chembl_id 속성이 붙어 있다
  flex_*_lig.pdb            유연 측쇄 도킹의 최적 포즈 (리간드)
  flex_belzutifan_sc.pdb    같은 도킹에서 움직인 측쇄 3개
  joi_bio.pml               장면 3개가 든 PyMOL 스크립트

무엇을 보면 되나
  M252 측쇄 뒤집힘이 이 묶음의 핵심이다. Ren et al. 2022 (Mol Pharmacol
  102:240) 가 belzutifan 결합 시 M252 가 포켓 밖으로 뒤집히고 H293 이
  수소결합, F254 가 pi-pi 상호작용을 한다고 보고했다. M252A 변이에서
  Ki 가 23 nM -> 3.5 uM 로 약 150배 악화된다.

  우리는 이 세 잔기를 유연하게 풀고 도킹했는데, 실제로 움직인 것은
  M252 하나였다. 측쇄 RMSD 2.37 A, CE 원자 3.96 A, SD 2.43 A 이동.
  리간드와의 접촉 거리가 1.53 A 에서 3.21 A 로 벌어진다.
  PHE254 는 0.16 A, HIS293 은 0.40 A 로 사실상 제자리다.

주의
  Vina 스코어는 순위용이다. 결합 자유에너지 추정에 쓰지 않는다.
  이 포즈들은 실험으로 확인된 바 없다.
  EPAS1 치료군은 BIAS_LIKELY 판정 — 상위권 순위의 75%가 구조 선택에서 나왔다.
  자세한 내용은 docs/context-of-use.md 참조.
