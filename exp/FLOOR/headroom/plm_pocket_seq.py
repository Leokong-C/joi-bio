"""포켓 PDB 285개 -> 서열 CSV. 중앙값 152잔기, 미지잔기 0."""
import os, glob, pandas as pd
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}
rows=[]
for p in sorted(glob.glob("casf/prot_trim/*.pdb")):
    seen=set(); seq=[]
    for ln in open(p):
        if ln.startswith("ATOM") and ln[12:16].strip()=="CA":
            k=(ln[21], ln[22:27])
            if k in seen: continue
            seen.add(k); seq.append(AA3.get(ln[17:20].strip(),"X"))
    rows.append({"code":os.path.basename(p)[:4],"n_res":len(seq),"seq":"".join(seq)})
pd.DataFrame(rows).to_csv("casf_pocket_seq.csv", index=False)
print(len(rows), "pockets")
