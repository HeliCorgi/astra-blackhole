"""Continuum threshold/domain audit for the KS positive-frequency sector.

Model:
    H^2 = -h^2 d_x^2 + exp(-2x),  H >= 0, on L2(R,dx).

The exact Liouville/Kontorovich-Lebedev generalized eigenfunctions are

    phi_k(x) = sqrt(2 k sinh(pi k))/pi * K_{ik}(exp(-x)/h), k>0,

with H phi_k = h k phi_k.  Their threshold expansion is

    phi_k(x) = sqrt(2/pi) k K_0(exp(-x)/h) + O(k^3),

pointwise for fixed x.  Therefore a state with threshold coefficient
B0=sqrt(2/pi) int K0(exp(-x)/h) psi(x) dx != 0 has spectral amplitude
c(k)=B0 k+O(k^3), so H^{-s}psi is locally square-integrable at k=0 iff
s<3/2.

The positive-frequency KG space is defined as the completion of D(H^1/2)
under ||Psi||_KG^2=(2/h)||H^1/2 Psi||^2.  The map
S=sqrt(h/2) H^(-1/2) extends by completion to a unitary map from L2 onto
that KG Hilbert space, even though H^(-1/2) is not a bounded L2 operator.

For the mass candidate, V=(i/h)[H,X] obeys -I <= V <= I in the continuum.
A convenient proof uses H V + V H = 2P and the Sylvester/semigroup formula
V=2 int_0^infty exp(-tH) P exp(-tH) dt together with |P|<=H.
This establishes I+V>=0 without a finite-box spectral gap.  It does NOT,
by itself, prove closability of
q_M[psi]=1/4 <H e^{X/2}psi,(I+V)H e^{X/2}psi>; the gap closes toward the
null/asymptotic channel, so that question remains separate.

Exploratory/model-internal mathematics only; not a singularity result.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import eig_banded
from scipy.special import k0


def continuum_packet(x: np.ndarray, h: float, *, x0: float = 2.0, p0: float = -1.0,
                     sigma: float | None = None, chirp: float = 0.0) -> np.ndarray:
    if sigma is None:
        sigma = math.sqrt(h)
    z = np.asarray(x, dtype=float) - x0
    norm = (2.0 * math.pi * sigma * sigma) ** (-0.25)
    return norm * np.exp(-z*z/(4.0*sigma*sigma) + 1j*(p0*z + chirp*z*z/2.0)/h)


def threshold_coefficient(*, h: float = 0.2, left: float = -4.0, right: float = 8.0,
                          dx: float = 0.005) -> complex:
    n = int(round((right-left)/dx))
    if n < 100 or not np.isclose(n*dx, right-left, atol=1e-12):
        raise ValueError("aligned threshold grid required")
    x = left + dx*np.arange(n+1)
    psi = continuum_packet(x, h)
    mode0 = math.sqrt(2.0/math.pi) * k0(np.exp(-x)/h)
    return complex(trapezoid(mode0*psi, x))


def threshold_power_for_H_minus_s(s: float) -> float:
    """Power p in int_0^eps k^p dk for |c(k)|~k and H^{-s}."""
    return 2.0 - 2.0*s


def threshold_domain_is_finite(s: float) -> bool:
    return threshold_power_for_H_minus_s(s) > -1.0


def finite_box_velocity_check(*, h: float = 0.2, left: float = -4.0,
                              right: float = 12.0, dx: float = 0.08) -> dict:
    """Independent finite-box check of the continuum commutator argument.

    Uses the repository's positive fourth-order Dirichlet discretization of H^2.
    P_comm=(i/2h)[H^2,X] is used so the Sylvester identity is tested without
    mixing in a separate finite-difference momentum stencil.
    """
    intervals = int(round((right-left)/dx))
    if intervals < 20 or not np.isclose(intervals*dx, right-left, atol=1e-12):
        raise ValueError("aligned box required")
    x = left + dx*np.arange(1, intervals)
    n = len(x)
    band = np.zeros((3, n))
    q = h*h/(12.0*dx*dx)
    band[0] = 30.0*q + np.exp(-2.0*x)
    band[0, [0, -1]] -= q
    band[1, :-1] = -16.0*q
    band[2, :-2] = q
    ev, U = eig_banded(band, lower=True, check_finite=True)
    H = (U*np.sqrt(ev)) @ U.T
    A = (U*ev) @ U.T
    X = np.diag(x)
    V = 1j/h*(H@X-X@H)
    V = (V+V.conj().T)/2.0
    P = 1j/(2.0*h)*(A@X-X@A)
    P = (P+P.conj().T)/2.0
    sylvester = np.linalg.norm(H@V+V@H-2.0*P)/max(np.linalg.norm(2.0*P), 1e-30)
    ve = np.linalg.eigvalsh(V)
    min_a_minus_p2 = float(np.linalg.eigvalsh(A-P@P)[0])
    return {
        "h": h,
        "left": left,
        "right": right,
        "dx": dx,
        "n": n,
        "sylvester_relative_residual": float(sylvester),
        "min_velocity_eigenvalue": float(ve[0]),
        "max_velocity_eigenvalue": float(ve[-1]),
        "min_I_plus_velocity_eigenvalue": float(1.0+ve[0]),
        "min_discrete_H2_minus_P2_eigenvalue": min_a_minus_p2,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    coeff_controls = [
        {"left": -4.0, "right": 8.0, "dx": 0.01},
        {"left": -6.0, "right": 10.0, "dx": 0.005},
        {"left": -8.0, "right": 12.0, "dx": 0.0025},
    ]
    coeff_rows = []
    for cfg in coeff_controls:
        b = threshold_coefficient(h=0.2, **cfg)
        coeff_rows.append({**cfg, "B0_real": float(b.real), "B0_imag": float(b.imag), "abs_B0": float(abs(b))})
    mags = [r["abs_B0"] for r in coeff_rows]
    rel_spread = (max(mags)-min(mags))/max(np.mean(mags), 1e-30)

    exponents = []
    for s in (0.5, 1.0, 1.5, 2.0):
        exponents.append({
            "s": s,
            "integrand_power": threshold_power_for_H_minus_s(s),
            "finite_at_zero_for_generic_nonzero_B0": threshold_domain_is_finite(s),
        })

    velocity_rows = [finite_box_velocity_check(right=r) for r in (8.0, 12.0, 16.0)]
    gaps = [r["min_I_plus_velocity_eigenvalue"] for r in velocity_rows]
    if not all(r["sylvester_relative_residual"] < 5e-11 for r in velocity_rows):
        raise RuntimeError("finite-box Sylvester identity regression")
    if not all(r["min_discrete_H2_minus_P2_eigenvalue"] > -1e-9 for r in velocity_rows):
        raise RuntimeError("discrete H^2-P^2 positivity regression")
    if not all(-1.0000001 <= r["min_velocity_eigenvalue"] <= r["max_velocity_eigenvalue"] <= 1.0000001
               for r in velocity_rows):
        raise RuntimeError("velocity escaped [-1,1]")
    if rel_spread > 1e-8 or min(mags) < 1e-6:
        raise RuntimeError(("threshold coefficient control failed", rel_spread, mags))

    result = {
        "schema": 1,
        "status": "PARTIAL",
        "model": "Vacuum KS positive-frequency Liouville Hamiltonian H=sqrt(-h^2 d_x^2+exp(-2x)) on L2(R,dx)",
        "units": "dimensionless repository conventions",
        "retained_variables": ["x", "positive-frequency T branch"],
        "omitted_variables": ["matter", "inhomogeneous modes", "environment"],
        "exact_continuum_spectral_representation": {
            "status": "PASS",
            "generalized_eigenfunction": "phi_k=sqrt(2 k sinh(pi k))/pi K_{i k}(exp(-x)/h), k>0",
            "energy": "E=h k",
            "spectrum": "continuous [0,infinity), no L2 zero mode",
            "threshold": "phi_k/k -> sqrt(2/pi) K_0(exp(-x)/h)",
        },
        "positive_frequency_KG_completion": {
            "status": "PASS",
            "space": "completion of D(H^(1/2)) in ||Psi||_KG^2=(2/h)||H^(1/2)Psi||_2^2",
            "same_state_map": "S=sqrt(h/2) H^(-1/2) extends by completion to a unitary L2 -> H_KG^+ map",
            "ordinary_L2_warning": "H^(-1/2) is unbounded on L2; not every KG amplitude produced by the completion is an ordinary L2 function.",
            "spectral_domain_rule": "D(H^a)={f: integral (h k)^(2a)|f_tilde(k)|^2 dk < infinity}; negative a uses the same spectral rule.",
        },
        "reference_packet_threshold": {
            "h": 0.2,
            "packet": "x0=2, p0=-1, sigma=sqrt(h), chirp=0",
            "B0_controls": coeff_rows,
            "relative_abs_B0_spread": float(rel_spread),
            "B0_nonzero": True,
            "spectral_amplitude": "c(k)=B0 k+O(k^3) at k->0",
            "H_minus_s_threshold_tests": exponents,
            "decision": "The reference packet is in ordinary-L2 D(H^-1/2) and D(H^-1), while a generic nonzero-B0 state reaches a logarithmic threshold obstruction at H^-3/2.",
        },
        "continuum_velocity": {
            "status": "PASS",
            "definition": "V=(i/h)[H,X] on a common core, equivalently the bounded Sylvester solution H V+V H=2P",
            "semigroup_formula": "V=2 int_0^infinity exp(-tH) P exp(-tH) dt",
            "bound": "-I <= V <= I because |P|<=H and ker(H)=0",
            "consequence": "I+V is a bounded positive operator/form in the continuum; finite-box positivity is not the only evidence.",
            "finite_box_controls": velocity_rows,
            "gap_decreases_with_box": bool(all(b<a for a,b in zip(gaps, gaps[1:]))),
            "coercive_gap_decision": "NO_UNIFORM_POSITIVE_GAP_ESTABLISHED",
        },
        "mass_continuum_form": {
            "status": "PENDING",
            "candidate": "q_M[psi]=(1/4)<H exp(X/2)psi,(I+V)H exp(X/2)psi>",
            "dense_core": "C_c^infinity(R) gives finite q_M",
            "what_is_now_proven": "The middle factor I+V is continuum-positive and bounded by 2.",
            "what_is_not_proven": "Closability/closedness of the full weighted nonlocal form and the associated self-adjoint mass operator.",
            "why_finite_box_gap_is_insufficient": "min spec(I+V) decreases toward zero under box expansion, consistent with the asymptotic null/left-moving channel; no coercive lower bound is available to transfer closure from an unweighted graph norm.",
        },
        "Y_clock_implication": {
            "status": "PENDING",
            "finding": "The same null-channel degeneracy explains why inverse-p_X is unbounded in the continuum even though selected finite-box packet norms can remain finite.",
        },
        "limitations": [
            "The normalized Liouville spectral representation and threshold expansion control the positive-frequency KS model only.",
            "The finite matrices are numerical checks of the analytic commutator argument, not a proof of mass-form closability.",
            "No quantum Kretschmann operator or black-hole singularity conclusion is promoted by this audit.",
        ],
        "reproducible_command": "python research/ks_continuum_domain_audit.py --out artifacts/ks-continuum/summary.json",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
