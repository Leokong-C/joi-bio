# Joi_bio — MET 키나제 ATP 사이트 (5HTI)
#
# PyMOL 에서:  pymol met.pml
#
# 주의: MET 은 대조군이 없어 게이트 판정이 ABSTAIN 이다.
#       이 포즈들로 어떤 결론도 내리지 않는다.

reinitialize
bg_color white
set ray_opaque_background, 0

load receptor_5HTI.pdb, met
load poses_rigid.sdf,   poses

hide everything
show cartoon, met
color grey80, met
set cartoon_transparency, 0.4, met

# 포즈 — 승인약과 상위 히트를 색으로 가른다
show sticks, poses
color yellow, poses and elem C
util.cnc("poses")

# 승인약 4개만 따로
select approved, poses and (resn UNL)
show sticks, poses

orient poses
zoom poses, 6

set stick_radius, 0.14
set cartoon_side_chain_helper, on

scene atp_site, store, message="MET ATP 사이트 (5HTI, 1.66 A). 포즈 19개 — 승인약 4 + 상위 히트 15"

print ""
print "=== Joi_bio MET ==="
print "구조  5HTI chain A, 2,304 원자, 공결정 리간드 66L (중원자 43)"
print "유사도 0.7195 (foretinib 과 거의 같은 화학형)"
print ""
print "판정  ABSTAIN — 대조군이 없다"
print ""
print "상위 22개 안에 활성이 17개 들어왔지만 승인약은 0개다."
print "  foretinib     30 위 (1.391%)"
print "  cabozantinib  33 위 (1.531%)"
print "  savolitinib  395 위 (18.321%)"
print "  crizotinib   690 위 (32.004%)"
print ""
print "승인약을 밀어낸 것은 중원자 25~34 의 작고 효율 높은 화합물들이다."
print "Vina 원점수가 작은 분자를 위로 올리는 경향과 맞지만, 원시 점수가 없어"
print "검증하지 못한 가설이다."
print ""
