"""ESM-2 포켓 임베딩이 바닥을 올리는가 — 바닥 6개 + 차원축소 대조 6개."""
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, GroupKFold, KFold
from scipy.stats import pearsonr
pm=pd.read_csv("casf2016_predictor_matrix.csv")[["code","logKa","target"]]
seq=pd.read_csv("casf_pocket_seq.csv"); E=np.load("casf_esm2_650M.npy")
emb=pd.DataFrame(E,columns=[f"esm{i}" for i in range(E.shape[1])]); emb.insert(0,"code",seq.code.values)
pk=pd.read_csv("headroom_features_pocket.csv"); lg=pd.read_csv("headroom_features_ligand.csv")
d=pm.merge(emb,on="code").merge(pk,on="code").merge(lg,on="code"); assert len(d)==285
y=d.logKa.values; g=d.target.values
ESM=[c for c in d if c.startswith("esm")]; PK=[c for c in pk.columns if c!="code"]; LG=[c for c in lg.columns if c!="code"]
def ev(mk,X,name,rows):
    o={"model":name}
    for s,cv,kw in [("disjoint",GroupKFold(10),{"groups":g}),("exposed",KFold(10,shuffle=True,random_state=0),{})]:
        o[s]=round(float(pearsonr(cross_val_predict(mk(),X,y,cv=cv,**kw),y)[0]),3)
    rows.append(o); print(f"  {name:34} disjoint {o['disjoint']:.3f}  exposed {o['exposed']:.3f}", flush=True)
RF=lambda: RandomForestRegressor(500,random_state=0,n_jobs=16)
X=lambda cols: np.nan_to_num(d[cols].to_numpy(np.float64),posinf=0,neginf=0)
r1=[]
for cols,name in [(PK,"ProteinBias (AA counts 21)"),(ESM,"ESM2-650M pocket (1280)"),(ESM+PK,"ESM2 + AA counts"),
                  (LG,"LigandBias (205)"),(LG+PK,"BothBias = Lig + AA counts"),(LG+ESM,"Lig + ESM2")]:
    ev(RF,X(cols),name,r1)
pd.DataFrame(r1).to_csv("headroom_esm2_floor.csv",index=False)
r2=[]; Xe=X(ESM); Xp=X(PK)
for k in [8,16,32,64]:
    ev(lambda k=k: make_pipeline(StandardScaler(),PCA(k,random_state=0),RF()),Xe,f"ESM2 -> PCA{k} -> RF",r2)
RG=lambda: make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,4,25)))
ev(RG,Xe,"ESM2 -> Ridge (linear)",r2); ev(RG,Xp,"AA counts -> Ridge (linear)",r2)
pd.DataFrame(r2).to_csv("headroom_esm2_control.csv",index=False)
