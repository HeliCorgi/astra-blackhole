"""Independent symbolic curvature calculation from the diagonal metric."""
from __future__ import annotations
from pathlib import Path
import json
import sympy as s

def main():
    t,x,th,ph=s.symbols('t x theta phi',real=True)
    a=s.Function('a')(t);b=s.Function('b')(t)
    coords=[t,x,th,ph];diag=[-s.S.One,a*a,b*b,b*b*s.sin(th)**2]
    G={}
    for i in range(4):
      for j in range(4):
       for k in range(4):
        v=(s.diff(diag[i],coords[k]) if i==j else 0)+(s.diff(diag[i],coords[j]) if i==k else 0)-(s.diff(diag[j],coords[i]) if j==k else 0)
        v=s.simplify(v/(2*diag[i]))
        if v!=0:G[i,j,k]=v
    gamma=lambda i,j,k:G.get((i,j,k),s.S.Zero)
    R={}
    for i in range(4):
     for j in range(4):
      for k in range(4):
       for l in range(4):
        val=s.diff(gamma(i,l,j),coords[k])-s.diff(gamma(i,k,j),coords[l])+sum(gamma(i,k,m)*gamma(m,l,j)-gamma(i,l,m)*gamma(m,k,j) for m in range(4))
        val=s.simplify(val)
        if val!=0:R[i,j,k,l]=val
    Ric=[s.simplify(sum(R.get((i,j,i,j),0) for i in range(4))) for j in range(4)]
    scalar=s.simplify(sum(Ric[i]/diag[i] for i in range(4)))
    Einstein=[s.simplify(Ric[i]/diag[i]-scalar/2) for i in range(4)]
    ha=s.diff(a,t)/a; hb=s.diff(b,t)/b
    expected=[-(2*ha*hb+hb**2+1/b**2),-(2*s.diff(b,t,2)/b+hb**2+1/b**2),-(s.diff(a,t,2)/a+s.diff(b,t,2)/b+ha*hb)]
    for i in range(3):assert s.simplify(Einstein[i]-expected[i])==0
    assert s.simplify(Einstein[2]-Einstein[3])==0
    K=s.simplify(sum(diag[i]*v*v/(diag[j]*diag[k]*diag[l]) for (i,j,k,l),v in R.items()))
    Kexpected=4*((s.diff(a,t,2)/a)**2+2*(s.diff(b,t,2)/b)**2+2*(ha*hb)**2+(hb**2+1/b**2)**2)
    assert s.simplify(K-Kexpected)==0
    output={'einstein_components_verified':True,'kretschmann_contraction_verified':True,
            'definition':'R^a_bcd and signature (-,+,+,+); metric computed directly.',
            'G_mixed':[str(e) for e in Einstein], 'K':str(K)}
    out=Path(__file__).parent/'results';out.mkdir(exist_ok=True)
    (out/'symbolic_geometry.json').write_text(json.dumps(output,indent=2))
    print('Direct Christoffel/Riemann computation verifies all implemented geometric formulas.')
if __name__=='__main__':main()
