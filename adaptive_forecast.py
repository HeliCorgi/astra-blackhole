"""Predictor-only extraction from the archived LTB benchmark; bodies unchanged.
Disagreement tolerance is not a bound on the true forecast error.
"""
from __future__ import annotations
import numpy as np

def fit_polynomials(time: np.ndarray, K: np.ndarray) -> dict:
    """Only past samples enter this function. No reference/model is accepted."""
    time = np.asarray(time, dtype=float)
    K = np.asarray(K, dtype=float)
    if time.ndim != 1 or time.shape != K.shape or len(time) < 8:
        raise ValueError('Need matching one-dimensional arrays with >=8 samples.')
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(K)):
        raise ValueError('Samples must be finite.')
    if not np.all(np.diff(time) > 0) or np.any(K <= 0):
        raise ValueError('Times must increase; curvatures must be positive.')
    width = float(time[-1]-time[0])
    x = (time-time[-1])/width
    q = (K/K[-1])**(-.25)
    coefficients = {}
    for degree in (1, 2, 3):
        design = np.column_stack([x**p for p in range(1, degree+1)])
        c = np.linalg.lstsq(design, q-1., rcond=None)[0]
        coefficients[degree] = np.r_[1., c]
    if coefficients[1][1] >= 0:
        raise ValueError('The linear fit does not predict growing curvature.')
    return {'time': float(time[-1]), 'width': width, 'K0': float(K[-1]),
            'coefficients': coefficients,
            'remaining_estimate': -width/coefficients[1][1]}

def ratio(fit: dict, degree: int, h: np.ndarray | float) -> np.ndarray:
    q = np.polynomial.polynomial.polyval(np.asarray(h)/fit['width'],
                                         fit['coefficients'][degree])
    if np.any(q <= 0):
        raise ValueError('Forecast polynomial reaches its spurious pole.')
    return q**(-4)

def predict(time: np.ndarray, K: np.ndarray, tolerance: float = .01) -> dict:
    """Quadratic forecast, shortening its horizon when orders 1/2/3 disagree.

    'tolerance' is an ensemble-disagreement threshold, NOT a proven error bound.
    The nominal horizon is half the linearly estimated remaining time.
    """
    fit = fit_polynomials(time, K)
    nominal_h = .5*fit['remaining_estimate']
    h = nominal_h
    for halves in range(20):
        try:
            hh = np.linspace(0., h, 17)
            estimates = np.array([ratio(fit, degree, hh) for degree in (1, 2, 3)])
            disagreement = float(np.max(estimates.max(axis=0)/estimates.min(axis=0)-1))
        except ValueError:
            disagreement = float('inf')
        if disagreement <= tolerance:
            return {**fit, 'nominal_h': nominal_h, 'adaptive_h': h,
                    'halvings': halves, 'disagreement': disagreement,
                    'linear_ratio': float(ratio(fit, 1, nominal_h)),
                    'quadratic_ratio': float(ratio(fit, 2, nominal_h)),
                    'adaptive_ratio': float(ratio(fit, 2, h))}
        h *= .5
    raise RuntimeError('Unable to choose a finite prediction horizon.')
