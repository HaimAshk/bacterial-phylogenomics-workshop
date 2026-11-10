import numpy as np, re, json
from scipy.stats import fisher_exact
K='/tmp/kit/'
# ---- tree
txt=open(K+'toy2_tree.nwk').read().strip().rstrip(';')
nodes=[]  # dict: name,parent,len,children
def parse(s):
    pos=0
    def rec(parent):
        nonlocal pos
        idx=len(nodes); nodes.append({'name':None,'parent':parent,'len':0.0,'ch':[]})
        if s[pos]=='(':
            pos+=1
            while True:
                c=rec(idx); nodes[idx]['ch'].append(c)
                if s[pos]==',': pos+=1; continue
                if s[pos]==')': pos+=1; break
        m=re.match(r'([A-Za-z0-9_.]*)(?::([0-9.eE-]+))?',s[pos:]); pos+=m.end()
        nodes[idx]['name']=m.group(1) or None
        nodes[idx]['len']=float(m.group(2)) if m.group(2) else 0.0
        return idx
    rec(None)
parse(txt)
N=len(nodes)
leaves=[i for i in range(N) if not nodes[i]['ch']]
# preorder
pre=[];st=[0]
while st:
    i=st.pop();pre.append(i);st.extend(reversed(nodes[i]['ch']))
post=pre[::-1]
# ---- data
rows={};cur=None
for l in open(K+'toy2_phyletic_pattern.fas'):
    l=l.strip()
    if l.startswith('>'):cur=l[1:];rows[cur]=''
    else:rows[cur]+=l
names=[nodes[i]['name'] for i in leaves]
G=np.array([[int(c) for c in rows[n]] for n in names])  # genomes x genes
truth=[l.rstrip('\n').split('\t') for l in open(K+'toy2_truth.tsv')][1:]
role=[t[2] for t in truth]; ann=[l.rstrip('\n').split('\t')[2] for l in open(K+'toy2_annotation.tsv')][1:]
# ---- trait: clade A pathogenic by descent; elsewhere pathogen iff carries M1
m1=[i for i,r in enumerate(role) if r=='M1']
carrier=G[:,m1].sum(1)>=len(m1)/2
trait=np.array([1 if (n.startswith('A') or carrier[k]) else 0 for k,n in enumerate(names)])
# ---- Fitch
leafidx={l:k for k,l in enumerate(leaves)}
def fitch(tips):
    S=[None]*N
    for i in post:
        if not nodes[i]['ch']: S[i]={int(tips[leafidx[i]])}
        else:
            a=[S[c] for c in nodes[i]['ch']]
            inter=set.intersection(*a)
            S[i]=inter if inter else set.union(*a)
    st=[None]*N
    for i in pre:
        if nodes[i]['parent'] is None: st[i]=min(S[i]) if len(S[i])==1 else 0
        else:
            p=st[nodes[i]['parent']]
            st[i]=p if p in S[i] else min(S[i])
    return np.array(st)
br=[i for i in range(N) if nodes[i]['parent'] is not None]
par=np.array([nodes[i]['parent'] if nodes[i]['parent'] is not None else 0 for i in range(N)])
blen=np.array([nodes[i]['len'] for i in range(N)])
def changes(st): 
    d=st-st[par]; d[0]=0; return d
def scores(tips,st_g,st_p):
    dg=changes(st_g);dp=changes(st_p)
    sim=abs(int((dg*dp)[br].sum()))
    same=(st_g==st_p)[br]
    sub=abs(2*(blen[br]*same).sum()/blen[br].sum()-1)
    g=tips-tips.mean();p=trait-trait.mean()
    term=abs((g*p).sum()/np.sqrt((g**2).sum()*(p**2).sum())) if g.any() else 0
    return term,sim,sub
stp=fitch(trait)
rng=np.random.default_rng(7)
def simulate(k,root):
    # choose k branches weighted by length; toggle downward
    w=blen[br]/blen[br].sum()
    ch=rng.choice(br,size=min(k,len(br)),replace=False,p=w)
    flip=np.zeros(N,int);flip[ch]=1
    st=np.zeros(N,int);st[0]=root
    for i in pre[1:]: st[i]=st[par[i]]^flip[i]
    return st
obs=[];ks=[];sts=[]
for j in range(G.shape[1]):
    st=fitch(G[:,j]);sts.append(st)
    ks.append(int(np.abs(changes(st)[br]).sum()));obs.append(scores(G[:,j],st,stp))
obs=np.array(obs)
NN=400;null=[]
for j in range(G.shape[1]):
    for _ in range(NN):
        s=simulate(max(ks[j],1),sts[j][0]); tips=s[leaves]
        if tips.min()==tips.max(): continue
        null.append(scores(tips,fitch(tips),stp))
null=np.array(null)
pv=np.array([[ (1+(null[:,t]>=obs[j,t]).sum())/(1+len(null)) for t in range(3)] for j in range(G.shape[1])])
fis=[]
for j in range(G.shape[1]):
    a=int(((G[:,j]==1)&(trait==1)).sum());b=int(((G[:,j]==1)&(trait==0)).sum());c=int(((G[:,j]==0)&(trait==1)).sum());d=int(((G[:,j]==0)&(trait==0)).sum())
    fis.append(fisher_exact([[a,b],[c,d]])[1])
fis=np.array(fis)
print('trait pathogens',trait.sum(),'carriers',carrier.sum())
print('naive Fisher p<0.01:', [(j,role[j]) for j in range(60) if fis[j]<0.01])
print('treeWAS simult p<0.01:', [(j,role[j],round(pv[j,1],4)) for j in range(60) if pv[j,1]<0.01])
print('any-score min p<0.01:', [(j,role[j]) for j in range(60) if pv[j].min()<0.01])
from collections import Counter
print(Counter(role))
for j in [0,1,2]: print(j,role[j],ks[j],obs[j],pv[j],fis[j])
np.savez('/home/claude/gwas/sim.npz',G=G,trait=trait,obs=obs,null=null,pv=pv,fis=fis,ks=np.array(ks),stp=stp,sts=np.array(sts),names=np.array(names))
json.dump({'role':role,'ann':ann},open('/home/claude/gwas/meta.json','w'))
