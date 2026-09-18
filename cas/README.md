# CAS gate contract

This directory defines the evidence contract for independent symbolic checks.
No xAct or Cadabra result is committed yet, so the GR reduction gate remains PENDING.

## Intended backends

- xAct/xTensor/xCoba for GR tensor and component calculations.
- Cadabra for an independent field-theory/tensor-algebra derivation.

These tools are not clean-clone dependencies of the repository at this stage.

## Required result manifest

A future CAS run must commit a machine-readable manifest containing at least:
backend/version, source/input and output SHA-256, metric/signature convention,
ADM/extrinsic-curvature convention, retained/omitted variables, lapse/shift,
boundary terms, reduced Lagrangian, canonical momenta, reduced Hamiltonian,
Hamiltonian constraint, Poisson-bracket convention, coefficient table,
implementation comparison, and explicit mismatches.

A configured script is not a passed gate. The manifest must point to an actually executed result.

## Promotion rule

The GR reduction gate may become PASS only after an independent CAS derivation has
been executed and its coefficients/signs/constraint structure have been compared
against the numerical model. A second backend or independent analytic certificate
should be used for high-confidence publication claims.

CAS success does not establish a unique physical quantization, clock, inner product,
or factor ordering.
