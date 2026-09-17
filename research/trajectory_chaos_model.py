"""Trajectory-vs-chaos controls for the reduced black-hole interior studies.

This file does not claim that a physical black-hole interior is chaotic or random.
It compares the exact one-degree-of-freedom Kantowski--Sachs principal-symbol
flow used by the reduced WDW study with the standard Gauss/BKL era map as a
known deterministic-chaos control.
"""
from __future__ import annotations
import numpy as np


def ks_energy(x, p):
    x, p = np.asarray(x, float), np.asarray(p, float)
    return np.sqrt(p*p + np.exp(-2*x))


def ks_path(x0, p0, t):
    """Exact flow of H=sqrt(p^2+exp(-2x))."""
    x0, p0, t = np.asarray(x0, float), np.asarray(p0, float), np.asarray(t, float)
    E = ks_energy(x0, p0)
    u0 = np.arcsinh(p0*np.exp(x0))
    u = u0 + t
    x = np.logaddexp(u, -u) - np.log(2.0) - np.log(E)
    p = E*np.tanh(u)
    return x, p


def finite_time_lyapunov_ks(x0=2.0, p0=-1.0, perturb=(0.0, 1e-8), times=None):
    """Euclidean phase-space finite-time separation exponent.

    It is coordinate/norm dependent at finite time; the long-time limit is the
    relevant chaos diagnostic. This is a numerical check, not a general proof.
    """
    if times is None:
        times = np.geomspace(1.0, 1e4, 81)
    times = np.asarray(times, float)
    dx0, dp0 = map(float, perturb)
    d0 = np.hypot(dx0, dp0)
    if d0 <= 0 or np.any(times <= 0):
        raise ValueError('Nonzero perturbation and positive times required.')
    x1, p1 = ks_path(x0, p0, times)
    x2, p2 = ks_path(x0+dx0, p0+dp0, times)
    d = np.hypot(x2-x1, p2-p1)
    return np.log(d/d0)/times, d


def gauss_map(x):
    """Gauss map, the standard era-map representative used in BKL statistics."""
    x = np.asarray(x, float)
    if np.any((x <= 0) | (x >= 1)):
        raise ValueError('Gauss-map state must lie in (0,1).')
    y = 1.0/x
    return y - np.floor(y)


def gauss_orbit(x0, steps):
    if steps < 1:
        raise ValueError('Positive number of steps required.')
    x = float(x0)
    out = np.empty(steps+1)
    out[0] = x
    for i in range(steps):
        x = float(gauss_map(x))
        if x == 0.0:
            raise ValueError('Rational/preperiodic orbit hit zero; use an irrational control.')
        out[i+1] = x
    return out


def gauss_finite_lyapunov(x0, steps):
    orbit = gauss_orbit(x0, steps)
    return float(np.mean(-2*np.log(orbit[:-1])))


def sample_gauss_invariant(n, seed=20260917):
    """Sample invariant density rho(x)=1/[ln2(1+x)] on (0,1)."""
    if n < 1:
        raise ValueError('Positive sample size required.')
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    return np.exp2(u)-1.0


def gauss_ensemble_lyapunov(n=200000, seed=20260917):
    x = sample_gauss_invariant(n, seed)
    return float(np.mean(-2*np.log(x)))


def gauss_exact_lyapunov():
    return float(np.pi*np.pi/(6*np.log(2.0)))
