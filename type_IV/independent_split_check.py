# independent_split_check.py -- independent check (does not import or copy biq.py / split_checks.py).
# 1) CM example (Kubota rank): a primitive CM type of M = L*C3 (C3 cyclic cubic, totally real) with
#    multiplicities (2,1) at sigma1 and (1,2) at sigma2; prints dim MT (Kubota rank), dim Hg = rank - 1 and
#    dim ker(U_M -> U_K) = 5 (dimension condition for almost-neatness of a member with End^0 = M larger than L).
# 2) split lemma: det formula, Hilbert norm test,
# hyperbolicity of 12-dim rational trace form via signature+discriminant+Hasse invariants; controls a=1, random a>>0.
import os, random, sys, json, itertools, math, time
from fractions import Fraction as Fr
import numpy as np
from cypari import pari
T0=time.time(); BUDGET=420
# ---------- CM check: M = L*C3 (C3 cyclic cubic totally real), Gal = {1,g,c,cg} x Z/3 ; bits (gbit,cbit)
G=[(h,k) for h in range(4) for k in range(3)]
mul=lambda x,y:(x[0]^y[0],(x[1]+y[1])%3)
Phi={(0,0),(0,1),(1,0),(2,2),(3,1),(3,2)}   # h: 0=sigma1,1=sigma2(g),2=conj sigma1(c),3=conj sigma2(cg)
cPhi={mul((2,0),x) for x in Phi}
assert not (Phi & cPhi) and (Phi|cPhi)==set(G)
rows=[[1 if mul(gg,x) in Phi else 0 for x in G] for gg in G]
rank=int(np.linalg.matrix_rank(np.array(rows)))
stab=[gg for gg in G if {mul(gg,x) for x in Phi}==Phi]
cm=dict(n_tau=sum(1 for x in Phi if x[0] in (0,1)), n_s1=sum(1 for x in Phi if x[0]==0), n_s2=sum(1 for x in Phi if x[0]==1),
        kubota_rank_dimMT=rank, dimHg=rank-1, dim_U_MK=5, stabilizer=stab)
print("CM check:",cm,flush=True)
# ---------- L arithmetic, basis 1,s,t,st ; s^2=m, t^2=-d
def Lmul(x,y,m,d):
    a1,b1,c1,e1=x; a2,b2,c2,e2=y
    return (a1*a2+m*b1*b2-d*c1*c2-m*d*e1*e2, a1*b2+b1*a2-d*(c1*e2+e1*c2),
            a1*c2+c1*a2+m*(b1*e2+e1*b2), a1*e2+e1*a2+b1*c2+c1*b2)
conj=lambda x:(x[0],x[1],-x[2],-x[3]); ONE=(1,0,0,0)
def psign(p):
    s=1
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            if p[i]>p[j]: s=-s
    return s
def Ldet(H,m,d):
    tot=[0,0,0,0]
    for p in itertools.permutations(range(3)):
        pr=ONE
        for i in range(3): pr=Lmul(pr,H[i][p[i]],m,d)
        s=psign(p)
        for k in range(4): tot[k]+=s*pr[k]
    return tuple(tot)
def emb(x,eps,m,d):
    r=math.sqrt(m)*eps; return complex(x[0]+x[1]*r, math.sqrt(d)*(x[2]+x[3]*r))
def sig(H,eps,m,d):
    ev=np.linalg.eigvalsh(np.array([[emb(H[i][j],eps,m,d) for j in range(3)] for i in range(3)]))
    if min(abs(ev))<1e-7: return None
    return (int(sum(ev>0)),int(sum(ev<0)))
def Fsign(p,q,m):
    if p>=0 and q>=0: return 0 if (p==0 and q==0) else 1
    if p<=0 and q<=0: return -1
    if p>0: return 1 if p*p>q*q*m else -1
    return 1 if q*q*m>p*p else -1
NF=lambda x,m:x[0]*x[0]-m*x[1]*x[1]
TrLK=lambda x:(2*x[0],2*x[2])
def Kmul(x,y,d): return (x[0]*y[0]-d*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def Kinv(x,d):
    n=x[0]*x[0]+d*x[1]*x[1]; return (Fr(x[0])/n, Fr(-x[1])/n)
def Kdet(M,d):
    M=[[(Fr(a),Fr(b)) for (a,b) in row] for row in M]; n=len(M); det=(Fr(1),Fr(0))
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]!=(0,0)),None)
        if p is None: return (Fr(0),Fr(0))
        if p!=c: M[c],M[p]=M[p],M[c]; det=(-det[0],-det[1])
        det=Kmul(det,M[c][c],d); inv=Kinv(M[c][c],d)
        for r in range(c+1,n):
            f=Kmul(M[r][c],inv,d)
            if f!=(0,0):
                newrow=[]
                for k in range(n):
                    q=Kmul(f,M[c][k],d); newrow.append((M[r][k][0]-q[0],M[r][k][1]-q[1]))
                M[r]=newrow
    return det
def hmat(H,a,m,d):
    U=[ONE,(0,1,0,0)]; M=[[None]*6 for _ in range(6)]
    for i in range(3):
        for u in range(2):
            for j in range(3):
                for v in range(2):
                    M[2*i+u][2*j+v]=TrLK(Lmul(Lmul(a,U[u],m,d),Lmul(U[v],H[i][j],m,d),m,d))
    return M
def gram12(M,d):   # B(x,y)=Tr_{K/Q} h(x,y), Q-basis w_k, t*w_k
    G=[[Fr(0)]*12 for _ in range(12)]
    for k in range(6):
        for l in range(6):
            x,y=M[k][l]; G[k][l]=Fr(2*x); G[6+k][l]=Fr(-2*d*y); G[k][6+l]=Fr(2*d*y); G[6+k][6+l]=Fr(2*d*x)
    assert all(G[i][j]==G[j][i] for i in range(12) for j in range(12))
    return G
def diag_form(G):
    n=len(G); G=[r[:] for r in G]; out=[]; idx=list(range(n))
    while idx:
        p=next((k for k in idx if G[k][k]!=0),None)
        if p is None:
            f=next(((k,l) for k in idx for l in idx if k!=l and G[k][l]!=0),None)
            if f is None: out+= [Fr(0)]*len(idx); break
            k,l=f
            for r in range(n): G[k][r]+=G[l][r]
            for r in range(n): G[r][k]+=G[r][l]
            p=k
        piv=G[p][p]; out.append(piv); rest=[k for k in idx if k!=p]
        for k in rest:
            f=G[k][p]/piv
            if f!=0:
                for r in rest: G[k][r]-=f*G[p][r]
        idx=rest
    return out
def primes_of(v):
    v=abs(int(v)); return [] if v<=1 else [int(p) for p in pari('factor(%d)[,1]'%v)]
def is_norm_K(D,d):
    D=Fr(D)
    if D<=0: return False
    n=D.numerator*D.denominator; ps={2}|set(primes_of(n))|set(primes_of(d))
    return all(int(pari.hilbert(n,-d,p))==1 for p in ps)
def hyperbolic12(diag):
    if any(x==0 for x in diag) or sum(1 for x in diag if x>0)!=6: return False
    ints=[x.numerator*x.denominator for x in diag]; prod=1
    for v in ints: prod*=v
    if prod<0 or math.isqrt(prod)**2!=prod: return False
    ps={2}
    for v in ints: ps|=set(primes_of(v))
    for p in ps:
        c=1
        for i in range(12):
            for j in range(i+1,12): c*=int(pari.hilbert(ints[i],ints[j],p))
        if c!=(-1 if p==2 else 1): return False
    return True
def test(H,a,m,d,Delta):
    M=hmat(H,a,m,d); det=Kdet(M,d)
    formula = det[1]==0 and det[0]==Fr((4*m)**3*NF(a,m)**3*NF(Delta,m))
    hil=is_norm_K(-det[0],d); witt=hyperbolic12(diag_form(gram12(M,d)))
    return formula,hil,witt
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 20261005)
ds=[1,2,3,5,6,7,11,15,19,23,31]; ms=[2,3,5,6,7,10,11,13,14,17,21]
st=dict(cases=0,a0_totpos=0,a0_formula=0,a0_hilbert=0,a0_witt=0,c1_split=0,c1_agree=0,c1_formula=0,cR_split=0,cR_agree=0,cR_formula=0,fails=[])
pairs=[(d,m) for d in ds for m in ms]; random.shuffle(pairs)
for rep in range(2):
  for (d,m) in pairs:
    if time.time()-T0>BUDGET: break
    r=lambda: random.randint(-3,3)
    while True:
        H=[[None]*3 for _ in range(3)]
        for i in range(3):
            H[i][i]=(r(),r(),0,0)
            for j in range(i+1,3):
                x=(r(),r(),r(),r()); H[i][j]=x; H[j][i]=conj(x)
        s1=sig(H,1,m,d); s2=sig(H,-1,m,d)
        if s1 and s2 and {s1,s2}=={(2,1),(1,2)}: break
    eps=1 if s1==(2,1) else -1
    Delta=Ldet(H,m,d); assert Delta[2]==0 and Delta[3]==0
    a0=Lmul((0,-eps,0,0),Delta,m,d)
    pos=Fsign(a0[0],a0[1],m)==1 and Fsign(a0[0],-a0[1],m)==1
    f,h,w=test(H,a0,m,d,Delta)
    y=random.randint(-6,6); aR=(int(abs(y)*math.sqrt(m))+random.randint(1,12),y,0,0)
    assert Fsign(aR[0],aR[1],m)==1 and Fsign(aR[0],-aR[1],m)==1
    f1,h1,w1=test(H,ONE,m,d,Delta); fR,hR,wR=test(H,aR,m,d,Delta)
    st['cases']+=1; st['a0_totpos']+=pos; st['a0_formula']+=f; st['a0_hilbert']+=h; st['a0_witt']+=w
    st['c1_split']+=h1; st['c1_agree']+=(h1==w1); st['c1_formula']+=f1
    st['cR_split']+=hR; st['cR_agree']+=(hR==wR); st['cR_formula']+=fR
    if not (pos and f and h and w and h1==w1 and hR==wR and f1 and fR): st['fails'].append(dict(d=d,m=m,H=str(H)))
print(json.dumps(st),"time %.0fs"%(time.time()-T0),flush=True)
json.dump(dict(cm=dict((k,str(v)) for k,v in cm.items()),split=st),open(sys.argv[2] if len(sys.argv)>2 else os.path.join(os.path.dirname(os.path.abspath(__file__)),'independent_split_check.json'),'w'))
