"""Algebraic concordance checks, not proofs of physical equivalence.

Read docs/LITERATURE_CONCORDANCE_ja.md for the source equations and scope.
No simulation, repository access, or network access is performed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as s


def run_checks() -> dict:
    checks = []

    def zero(name, expression):
        reduced = s.simplify(expression.rewrite(s.exp))
        if reduced != 0:
            raise AssertionError(f'{name}: residual={reduced}')
        checks.append({'name': name, 'passed': True, 'residual': str(reduced)})

    x, tau = s.symbols('x tau', real=True)
    h = s.Symbol('h', positive=True)
    Om, be = s.symbols('Omega beta', real=True)
    f = s.Function('f')
    om_map = (x + s.log(4))/s.sqrt(3)
    be_map = tau/s.sqrt(3)
    F = f(s.sqrt(3)*Om-s.log(4), s.sqrt(3)*be)
    old = -h*h*s.diff(F,Om,2) + h*h*s.diff(F,be,2) + 48*s.exp(-2*s.sqrt(3)*Om)*F
    target = h*h*(s.diff(f(x,tau),tau,2)-s.diff(f(x,tau),x,2))+s.exp(-2*x)*f(x,tau)
    zero('LDO2009_Eq4_to_repo_WDW', old.subs({Om:om_map,be:be_map}).doit()/3-target)

    U, V = s.symbols('U V', real=True)
    phi = s.Function('phi')(U,V)
    DT = lambda z: -U*s.diff(z,U) + V*s.diff(z,V)
    DX = lambda z: -U*s.diff(z,U) - V*s.diff(z,V)
    diff = DT(DT(phi))-DX(DX(phi))
    zero('Chiba_chain_rule', diff+4*U*V*s.diff(phi,U,V))
    um = -s.exp(-x-tau)/4
    vm = s.exp(-x+tau)/4
    zero('Chiba_mixed_derivative_coefficient',(-4*U*V).subs({U:um,V:vm})-s.exp(-2*x)/4)
    # Both operators act on the same scalar function; this does not construct
    # a unitary map between the chosen physical Hilbert spaces.
    kappa = h/2
    zero('Chiba_Eq14_parameter_normalization',1/kappa**2-4/h**2)
    # Negative control: kappa=h, instead of h/2, must not be accepted.
    wrong = s.simplify(1/h**2-4/h**2)
    if wrong == 0:
        raise AssertionError('Wrong normalization escaped detection')
    checks.append({'name':'wrong_kappa_negative_control','passed':True,'nonzero_residual':str(wrong)})
    zero('areal_radius_mapping', -um-s.exp(-x-tau)/4)
    zero('radial_scale_mapping', -vm/um-s.exp(2*tau))
    tm = (um+vm)/2
    xm = (vm-um)/2
    zero('minisuperspace_clock_is_boost_like',s.expand(tm-s.exp(-x)*s.sinh(tau)/4))
    jac = s.Matrix([[s.diff(tm,tau),s.diff(tm,x)], [s.diff(xm,tau),s.diff(xm,x)]])
    g = s.simplify(jac.T*s.diag(-1,1)*jac)
    if s.simplify(g - s.exp(-2*x)/16*s.diag(-1,1)) != s.zeros(2):
        raise AssertionError('Minisuperspace metric mismatch')
    checks.append({'name':'minisuperspace_metric_pullback','passed':True,'residual':'zero 2x2 matrix'})

    # A finite-mode demonstration of the SAME-clock KG-to-L2 rescaling.
    # This is not an infinite-dimensional domain theorem or a clock-change map.
    e1,e2,z1,z2=s.symbols('e1 e2 z1 z2',positive=True)
    kg=2*e1*z1*z1+2*e2*z2*z2
    mapped=(s.sqrt(2*e1)*z1)**2+(s.sqrt(2*e2)*z2)**2
    zero('same_clock_KG_to_L2_finite_modes',kg-mapped)
    if s.simplify(kg-(z1*z1+z2*z2)) == 0:
        raise AssertionError('Untransformed wavefunction norms were equated')
    checks.append({'name':'same_wavefunction_norms_not_identical','passed':True})

    # CHS2023 Eq29, classical massless axisymmetric specialization:
    # f(q,mu)=exp[-q*A(mu)] => integral q^3 f dq = 6/A(mu)^4.
    q,A=s.symbols('q A',positive=True)
    zero('CHS2023_energy_integrated_entropy_ansatz',s.integrate(q**3*s.exp(-A*q),(q,0,s.oo))-6/A**4)
    y,a,b=s.symbols('y a b',real=True)
    second=s.simplify(s.diff(-4*s.log(a+b*y),y,2))
    zero('thermal_ME_log_curvature',second-4*b*b/(a+b*y)**2)
    if second == 0:
        raise AssertionError('Different entropy closures were equated')
    checks.append({'name':'angular_exponential_not_thermal_ME','passed':True,'difference':'d2(log I)/d(mu^2)^2 = 4*b^2/(a+b*mu^2)^2, not zero for b != 0'})

    # Only the elementary Abel integral step, NOT the spectral remainder.
    p,X,T=s.symbols('p X T',positive=True)
    elementary=2*p*X/(X**2+p**2)**2
    leading=s.limit(X**3*elementary.subs(p,s.I*T),X,s.oo)
    zero('conditional_threshold_integral_leading_term',leading-2*s.I*T)
    return {'schema':1, 'comparison_commit':'2d309769ca155617f95cd399b0844a382c4556d0',
            'sympy_version':s.__version__,'passed_checks':len(checks),'checks':checks,
            'scope':'Finite symbolic identities and negative controls only. No Lean proof, PDE simulation, spectral remainder estimate, Hilbert-space equivalence or empirical test.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    text=json.dumps(run_checks(),ensure_ascii=False,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(text,encoding='utf-8')
    print(text,end='')

if __name__=='__main__':
    main()
