"""Independent NumPy/SciPy checks of the algebra, not execution of MATLAB code."""
import itertools, json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular
from scipy.fft import ifftn, ifftshift, fftshift

rng=np.random.default_rng(2026)
records=[]
def record(name,actual,expected,tolerance):
    error=float(np.max(np.abs(actual-expected)))
    assert error <= tolerance,(name,error,tolerance)
    records.append(dict(check=name,max_absolute_error=error,tolerance=tolerance))
def pad(K,target):
    delta=np.array(target)-K.shape[:3]
    out=np.zeros(tuple(target)+K.shape[3:],dtype=K.dtype)
    out[tuple(slice(int(d//2),int(d//2+n)) for d,n in zip(delta,K.shape[:3]))]=K
    return out
def fft(K,axes):
    return fftshift(ifftn(ifftshift(K,axes=axes),axes=axes),axes=axes)*np.sqrt(np.prod([K.shape[a] for a in axes]))

shape=(9,7,3,2,3,4)
K=(rng.standard_normal(shape)+1j*rng.standard_normal(shape)).astype(np.complex64)
noise=(rng.standard_normal((4,400))+1j*rng.standard_normal((4,400))).astype(np.complex64)
L=np.linalg.cholesky(noise@noise.conj().T/399+np.eye(4,dtype=np.float32)*1e-6)
X=K.reshape(-1,4,order='F').T
full=solve_triangular(L,X,lower=True)
chunk=np.empty_like(full[[3,1]])
for first in range(0,X.shape[1],17):
    chunk[:,first:first+17]=solve_triangular(L,X[:,first:first+17],lower=True)[[3,1]]
record('full-channel whitening then output selection',chunk,full[[3,1]],3e-6)
Kw=full[[3,1]].T.reshape(shape[:5]+(2,),order='F')
for target,axes in [((14,11,3),(0,1)),((14,11,5),(0,1,2)),((9,7,3),(0,1,2))]:
    expected=fft(pad(Kw,target),axes)
    actual=np.empty_like(expected)
    for t,f in itertools.product(range(shape[3]),range(shape[4])):
        actual[:,:,:,t:t+1,f:f+1,:]=fft(pad(Kw[:,:,:,t:t+1,f:f+1,:],target),axes)
    record(f'blocked spatial FFT {target}/{axes}',actual,expected,1e-5)
    delta=np.array(target)-np.array(shape[:3])
    twice=np.pad(Kw,[(int(d//2),0) for d in delta]+[(0,0)]*3)
    twice=np.pad(twice,[(0,int(d-d//2)) for d in delta]+[(0,0)]*3)
    record(f'one-allocation centred padding {target}',pad(Kw,target),twice,0)

def llr(X,new):
    H,W,N=X.shape; p=5;step=3;thr=.1*.5*np.sqrt(max(p*p,N));Y=X.copy()
    wgt=np.zeros((H,W),np.float32) if new else np.zeros_like(X)
    if new:
        for i,j in itertools.product(range(0,H-p+1,step),range(0,W-p+1,step)):
            wgt[i:i+p,j:j+p]+=1
    for _ in range(2):
        acc=np.zeros_like(X)
        if not new: wgt[:]=0
        for i,j in itertools.product(range(0,H-p+1,step),range(0,W-p+1,step)):
            patch=Y[i:i+p,j:j+p,:].reshape(p*p,N,order='F')
            mean=patch.mean(0); U,S,V=np.linalg.svd(patch-mean,full_matrices=False)
            shrunk=np.maximum(S-thr,0)
            P=U@(shrunk[:,None]*V if new else np.diag(shrunk)@V)+mean
            acc[i:i+p,j:j+p,:]+=P.reshape(p,p,N,order='F')
            if not new:wgt[i:i+p,j:j+p,:]+=1
        if new:
            mask=wgt>0;Y[mask]=acc[mask]/wgt[mask,None]
        else:
            mask=wgt>0;Y[mask]=acc[mask]/wgt[mask]
        Y=np.maximum(Y,0)
    return Y
for frames in (7,72):
    source=np.abs(rng.standard_normal((11,13,frames))).astype(np.float32)
    record(f'LLR 2D coverage/diagonal-free shrinkage N={frames}',llr(source,True),llr(source,False),3e-6)

volume=np.arange(5*7*3*2).reshape(5,7,3,2)
for order in itertools.permutations(range(3)):
    for index in (0,volume.shape[order[2]]-1):
        permutation=order+(3,)
        expected=np.transpose(volume,permutation)[:,:,index,:]
        sel=[slice(None)]*4;sel[order[2]]=slice(index,index+1)
        actual=np.transpose(volume[tuple(sel)],permutation)[:,:,0,:]
        record(f'slice-before-permute {order}/{index}',actual,expected,0)

# Explicit field then time means and streaming reduction preserve axis semantics.
damaged=K.copy();damaged[0,1,2,0,1,1]=np.nan
reference=np.nanmean(np.nanmean(damaged,axis=4,keepdims=True),axis=3,keepdims=True)
total=np.zeros(shape[:3]+(1,1,4),np.complex64);nt=np.zeros(total.shape,np.float32)
for t in range(shape[3]):
    sums=np.zeros_like(total);nf=np.zeros(total.shape,np.float32)
    for f in range(shape[4]):
        block=damaged[:,:,:,t:t+1,f:f+1,:]; valid=~np.isnan(block)
        sums+=np.where(valid,block,0);nf+=valid
    means=sums/nf;valid=~np.isnan(means);total+=np.where(valid,means,0);nt+=valid
record('streaming field-then-time ACS average',total/nt,reference,3e-6)

report={'scope':'Independent algebra checks; MATLAB app and MATLAB tests were not executed.',
        'checks_passed':len(records),'checks':records}
out=Path(__file__).resolve().parent
out.mkdir(exist_ok=True)
(out/'independent_numeric_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'checks_passed':len(records),'max_error':max(r['max_absolute_error'] for r in records)}))
