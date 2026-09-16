"""Symbolic local derivatives of K from the moment hierarchy, no time fitting."""
from __future__ import annotations
from pathlib import Path
import sympy as s
import numpy as np
from scipy.special import beta

class CurvatureJets:
    def __init__(self, maximum: int=3, export: Path | None=None):
        ha,hb,B=s.symbols('Ha Hb B')
        self.ha,self.hb,self.B=ha,hb,B
        self.C={(i,k-i):s.Symbol(f'C{i}_{k-i}') for k in range(maximum+2) for i in range(k+1)}
        pr=self.C[1,0];pt=self.C[0,1]/2
        dhb=-(3*hb**2+B+pr)/2
        dha=-pt-dhb-ha**2-hb**2-ha*hb
        flow={ha:dha,hb:dhb,B:-2*hb*B}
        for (i,j),v in self.C.items():
            k=i+j
            if k<=maximum:
                flow[v]=-((1+2*i)*ha+(2+2*j)*hb)*v+(2*k-1)*(ha*self.C[i+1,j]+hb*self.C[i,j+1])
        K=s.expand(4*((dha+ha**2)**2+2*(dhb+hb**2)**2+2*(ha*hb)**2+(hb**2+B)**2))
        expr=[K]
        for _ in range(maximum):
            expr.append(s.expand(sum(s.diff(expr[-1],v)*flow[v] for v in expr[-1].free_symbols)))
        self.expressions=tuple(expr)  # Exact-algebra audits may substitute rational moments.
        self.symbols=[ha,hb,B]+list(self.C.values())
        self.functions=[s.lambdify(self.symbols,e,modules='numpy',cse=True) for e in expr]
        if export:
            export.write_text('\n\n'.join(f'd^{i}K/dt^{i} = {e}' for i,e in enumerate(expr)))

    def evaluate(self, basis, weights, geometry):
        moments=basis.moment_rows(5)[1:]@weights
        la,lb,ha,hb=geometry
        if abs(la)+abs(lb)>1e-12:raise ValueError('This evaluation assumes initially isotropic a=b=1.')
        vals=[ha,hb,1.]
        for i,j in self.C:
            vals.append(.5*beta(i+.5,j+1)*moments[i+j])
        return np.array([f(*vals) for f in self.functions],dtype=float)
