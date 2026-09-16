"""Positive ANGULAR-entropy closure, same evolved moments as closure_model.
Copyright 2026 HeliCorgi. Apache-2.0.

I/rho = exp(sum lambda_l P_l)/<exp(sum lambda_l P_l)>, l=2 or 2,4.
The multipliers solve an algebraic moment problem at each evaluation; they
are not additional dynamical variables or fitted physical parameters.
This maximizes angular Shannon entropy with measure dmu/2. It is NOT a
claim of maximizing the entropy of the full relativistic f(q,mu).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial import legendre as leg
from closure_model import ClosureModel
from collision_model import transport_matrix


class MomentDomainError(ValueError):
    """Input moments have no strictly positive interior realization."""


@dataclass
class AngularMaxEntropy:
    lmax: int
    quadrature_order: int = 80
    tolerance: float = 2e-13
    max_iterations: int = 80

    def __post_init__(self):
        if self.lmax not in (2, 4):
            raise ValueError('Retain c0,c2 or c0,c2,c4 only.')
        if not isinstance(self.quadrature_order, int) or self.quadrature_order < 16:
            raise ValueError('Quadrature order must be an integer >=16.')
        if not np.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError('Positive finite nonlinear tolerance required.')
        if not isinstance(self.max_iterations, int) or self.max_iterations < 1:
            raise ValueError('Positive iteration limit required.')
        self.mu, w = leg.leggauss(self.quadrature_order)
        self.weights = w/2
        self.degrees = np.arange(2, self.lmax+1, 2)
        self.features = leg.legvander(self.mu, self.lmax)[:, self.degrees]
        self.tail_feature = leg.legval(self.mu, [0.]*(self.lmax+2)+[1.])

    def _statistics(self, multipliers):
        exponent = self.features @ multipliers
        shift = float(exponent.max())
        unnormalized = self.weights*np.exp(exponent-shift)
        partition = float(unnormalized.sum())
        probabilities = unnormalized/partition
        mean = probabilities @ self.features
        centered = self.features-mean
        covariance = (centered.T*probabilities) @ centered
        return mean, covariance, shift+np.log(partition), probabilities

    def fit(self, coefficients: np.ndarray) -> dict:
        c = np.asarray(coefficients, dtype=float)
        if c.shape != (1+self.lmax//2,) or not np.all(np.isfinite(c)) or c[0] <= 0:
            raise MomentDomainError('Finite moments and rho>0 required.')
        z = c/c[0]
        xmean = 1/3+2*z[1]/15
        if not 0 < xmean < 1:
            raise MomentDomainError('c2 outside the interior positive-moment cone.')
        if self.lmax == 4:
            xsquare = 1/5+4*z[1]/35+8*z[2]/315
            if not xmean*xmean < xsquare < xmean:
                raise MomentDomainError('c4 outside the interior positive-moment cone.')
        target = z[1:]/(2*self.degrees+1)
        # No path-dependent cache: identical moments always give the same fit.
        multipliers = np.zeros(len(target))
        for iteration in range(self.max_iterations):
            mean, covariance, logZ, probabilities = self._statistics(multipliers)
            gradient = mean-target
            error = float(np.max(abs(gradient)))
            if error <= self.tolerance:
                return dict(multipliers=multipliers, logZ=float(logZ),
                    normalized_residual=error, iterations=iteration,
                    hessian_condition=float(np.linalg.cond(covariance)),
                    tail=float((2*(self.lmax+2)+1)*c[0]*(probabilities@self.tail_feature)))
            try:
                step = np.linalg.solve(covariance, gradient)
            except np.linalg.LinAlgError as exc:
                raise RuntimeError('Singular entropy Hessian; no projection/clipping performed.') from exc
            descent = float(gradient@step)
            objective = logZ-float(multipliers@target)
            for halving in range(50):
                factor = 2.**(-halving)
                trial = multipliers-factor*step
                nm, _, nlz, _ = self._statistics(trial)
                new_error = float(np.max(abs(nm-target)))
                # Residual safeguard prevents near-optimum roundoff in logZ.
                decrease = nlz-float(trial@target) <= objective-1e-4*factor*descent
                if new_error < error*(1-1e-4*factor) or (decrease and error>1e-7):
                    multipliers = trial
                    break
            else:
                raise RuntimeError(f'Moment solve stalled with residual {error}; not accepted.')
        raise RuntimeError('Moment solve iteration limit reached; not accepted.')

    def log_density(self, fit: dict, mu: np.ndarray) -> np.ndarray:
        mu = np.asarray(mu, dtype=float)
        if not np.all(np.isfinite(mu)) or np.any(abs(mu)>1):
            raise ValueError('Direction cosines must lie in [-1,1].')
        return leg.legvander(mu, self.lmax)[:, self.degrees]@fit['multipliers']-fit['logZ']

    def minimum_log_density(self, fit: dict) -> float:
        """All polynomial stationary points, not a grid-only positivity check."""
        poly = np.zeros(self.lmax+1)
        poly[self.degrees] = fit['multipliers']
        candidates = [-1., 1.]
        for r in leg.legroots(leg.legder(poly)):
            if abs(complex(r).imag)<1e-10 and -1<=float(np.real(r))<=1:
                candidates.append(float(np.real(r)))
        return float(min(leg.legval(candidates, poly))-fit['logZ'])

    def independent_residual(self, coefficients: np.ndarray, fit: dict,
                             quadrature_order: int = 160) -> dict:
        """A second quadrature checks the SAME lambda, not a refitted answer."""
        mu, w = leg.leggauss(quadrature_order)
        density = np.exp(self.log_density(fit, mu))
        all_degrees = np.arange(0, self.lmax+3, 2)
        recovered = ((w/2*density)@leg.legvander(mu, self.lmax+2)[:, all_degrees])*(2*all_degrees+1)
        target = np.asarray(coefficients)/coefficients[0]
        return dict(moment_residual=float(max(abs(recovered[:-1]-target))),
                    tail_residual=float(abs(recovered[-1]-fit['tail']/coefficients[0])))


class PositiveClosureModel(ClosureModel):
    """Only the FIRST OMITTED moment replaces zero in the same transport rows."""
    def __init__(self, lmax: int, rate: float, quadrature_order: int = 80,
                 tolerance: float = 2e-13):
        super().__init__(lmax, rate)
        self.reconstruction = AngularMaxEntropy(lmax, quadrature_order, tolerance)
        size = lmax//2+1
        self.tail_coupling = transport_matrix(8)[:size, size].copy()

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        original = super().rhs(t, y)
        fit = self.reconstruction.fit(y[5:])
        original[5:] += (y[2]-y[3])*self.tail_coupling*fit['tail']
        return original
