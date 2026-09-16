"""Positive angular-distribution counterexamples on the existing Vlasov solver.

Units c=8*pi*G=1. Initial a=b=1. This is classical, homogeneous,
collisionless matter; no new force, closure law, or quantum effect is added.
"""
from __future__ import annotations
import numpy as np
import sympy as sp
from scipy.special import eval_legendre
from kinetic_model import MomentumBasis


class AngularBasis(MomentumBasis):
    """Initial angular number density (1 + amplitude*P_ell(mu))/2.

    mu is the INITIAL physical radial direction cosine, transported via the
    two conserved momentum variables of MomentumBasis, not re-isotropized.
    Only even ell>=4 are admitted to keep reflection symmetry and initial T.
    """
    def __init__(self, centers, ell: int = 4, amplitude: float = .8, **kwargs):
        if not isinstance(ell, int) or ell < 4 or ell % 2:
            raise ValueError('ell must be an even integer >= 4.')
        if not np.isfinite(amplitude) or abs(amplitude) >= 1:
            raise ValueError('Finite |amplitude| < 1 ensures strict positivity.')
        super().__init__(centers, **kwargs)
        if self.angular_order < ell + 2:
            raise ValueError('Angular quadrature too small for matching checks.')
        self.ell, self.amplitude = ell, amplitude
        self.aw = self.aw * (1 + amplitude * eval_legendre(ell, self.mu))
        self.quadw = self.rw[:, None] * self.aw[None, :]


def radial_weights(basis, density=.8):
    """Same positive radial spectrum for both angular states, never fitted to futures."""
    w = 1 / np.sqrt(1 + basis.centers**2)
    return w * density / (basis.moment_rows(0)[1] @ w)


def general_jets(jets, basis, weights, geometry):
    """Evaluate existing symbolic K derivatives with actual angular moments.

    CurvatureJets.evaluate assumes isotropy; deliberately do NOT call it here.
    """
    la, lb, ha, hb = geometry
    order = max(i+j for i,j in jets.C)
    c = basis.tensors(geometry, weights, order)
    values = [ha, hb, np.exp(-2*lb)] + [c[key] for key in jets.C]
    return np.array([f(*values) for f in jets.functions], dtype=float)


def angular_identities(ell: int):
    """Exact rational angular integrals (computer algebra, NOT a Lean proof)."""
    mu = sp.Symbol('mu')
    p = sp.legendre(ell, mu)
    out = {}
    for k in range(ell//2 + 1):
        for i in range(k+1):
            j = k-i
            value = sp.integrate(mu**(2*i)*(1-mu**2)**j*p, (mu,-1,1))/2
            value = sp.factor(value)
            if 2*k < ell and value != 0:
                raise AssertionError('Angular orthogonality failed.')
            out[f'{i},{j}'] = str(value)
    return out
