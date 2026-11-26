"""Floating-point sanity (NOT a proof): sample H=(1+g)/(1-g) in P_2, g Schur with g(0)=0, check
sum_{n<=5} r^n ||H_n||_1 <= B = 0.9931560506 at r=189/625. Expects max found < B."""
import numpy as np, random, sys
r=189/625; B=0.9931560506; D=6
def mul(a,b):
    out=np.zeros((D,D),complex)
    for i in range(D):
        for j in range(D-i):
            if a[i,j]!=0:
                out[i:,j:]+=a[i,j]*b[:D-i,:D-j]
    out[np.add.outer(np.arange(D),np.arange(D))>=D]=0
    return out
def inv1(a):  # 1/(1-a), a(0)=0
    out=np.zeros((D,D),complex); out[0,0]=1; p=out.copy()
    for _ in range(D):
        p=mul(p,a); out+=p
    return out
def val(g):
    one=np.zeros((D,D),complex); one[0,0]=1
    num=one+g; den=one-g
    h=mul(num,inv1(g))   # (1+g)/(1-g)
    tot=0
    for i in range(D):
        for j in range(D):
            if 1<=i+j<=5: tot+=abs(h[i,j])*r**(i+j)
    return tot
def make(par):
    p,th1,th2,c_re,c_im,t,k=par
    z1=np.zeros((D,D),complex); z1[1,0]=1
    z2=np.zeros((D,D),complex); z2[0,1]=1
    s=p*np.exp(1j*th1)*z1+(1-p)*np.exp(1j*th2)*z2
    c=0.95*np.tanh(abs(complex(c_re,c_im)))*np.exp(1j*np.angle(complex(c_re,c_im)))
    one=np.zeros((D,D),complex); one[0,0]=1
    phi=mul(s,mul(s-c*one,inv1(np.conj(c)*s)))   # s(s-c)/(1-conj(c)s): self-map of disc, phi(0)=0
    phi=np.exp(1j*t)*phi
    return phi
best=0
rng=np.random.default_rng(1)
def f(x): return val(make(np.array([0.5+0.5*np.tanh(x[0]),*x[1:6],0])))
for trial in range(40):
    x=rng.normal(size=6); fx=f(x); step=0.5
    for it in range(400):
        y=x+step*rng.normal(size=6); fy=f(y)
        if fy>fx: x,fx=y,fy
        else: step=max(step*0.995,0.01)
    best=max(best,fx)
print('max found',best,'B',B,'ok' if best<B else 'VIOLATION')
