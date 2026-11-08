"""ESM-2 650M 평균풀링 임베딩. CPU 16스레드 191초. 가중치는 SHA-256 검증."""
import hashlib, os, sys, time, numpy as np, pandas as pd, torch
sys.path.insert(0, os.path.abspath("vendor"))   # PLM.md §1 의 wheel 우회
import esm
SHA={"models/esm2_t33_650M_UR50D.pt":"ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c",
     "regression/esm2_t33_650M_UR50D-contact-regression.pt":"8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575"}
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for c in iter(lambda: f.read(1<<20), b""): h.update(c)
    return h.hexdigest()
ck=os.path.join(torch.hub.get_dir(),"checkpoints"); os.makedirs(ck,exist_ok=True)
for rel,dig in SHA.items():
    dst=os.path.join(ck,os.path.basename(rel))
    if os.path.isfile(dst) and sha(dst)!=dig: os.rename(dst,dst+".rejected")
    if not os.path.isfile(dst):
        torch.hub.download_url_to_file("https://dl.fbaipublicfiles.com/fair-esm/"+rel,dst,hash_prefix=dig,progress=False)
torch.set_num_threads(16)
df=pd.read_csv("casf_pocket_seq.csv")
model,alphabet=esm.pretrained.esm2_t33_650M_UR50D(); model=model.eval(); bc=alphabet.get_batch_converter()
out=[]; t0=time.time()
with torch.no_grad():
    for _,r in df.iterrows():
        _,_,toks=bc([(r.code,r.seq)])
        out.append(model(toks,repr_layers=[33])["representations"][33][0,1:-1].mean(0).numpy())
np.save("casf_esm2_650M.npy", np.vstack(out).astype(np.float32))
print("done", len(out), f"{time.time()-t0:.0f}s")
