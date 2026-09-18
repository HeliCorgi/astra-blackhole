# CAS gate contract

This directory defines the evidence contract for independent symbolic checks.

Bianchi IX now has two independently executed backends: pinned Cadabra 2.5.14 and Ubuntu Maxima/ctensor 5.46.0. Their results are stored under `cas/results/` and documented in `docs/BIANCHI_IX_GR_REDUCTION_CAS_ja.md`. The xAct source remains committed but unexecuted because repository CI has no licensed Wolfram Engine/xAct runtime. Two executed independent backends are sufficient for the current GR-reduction obligation, so the aggregate Bianchi IX gate is PASS.

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


## Current Bianchi IX status

- Cadabra 2.5.14: PASS
  - independent Euler-angle SU(2) spatial metric
  - component Christoffel/Ricci scalar
  - ADM kinetic invariant
  - Legendre transform to the Hamiltonian constraint
  - comparison against classical and quantum repository implementations
- Maxima/ctensor 5.46.0: PASS
  - independent ctensor curvature calculation
  - independent ADM / canonical algebra
  - repository implementation comparison
- xAct/xCoba source: PENDING optional third execution
- aggregate GR_REDUCTION: PASS

The executed Cadabra path uses Cadabra's SymPy scalar backend inside the pinned
Cadabra process. This is recorded explicitly rather than described as Cadabra's
abstract tensor `evaluate` path.


## Backend independence note

Cadabra 2.5.14 is executed through `cadabra2-cli`; its final scalar-component
calculation uses Cadabra's SymPy scalar backend. The second backend therefore does
not use SymPy: Maxima 5.46.0 uses the Maxima `ctensor` package for the Euler-angle
metric curvature and Maxima's own algebra for the ADM and Legendre steps.

This distinction is why Maxima is counted as an independent second backend.
xAct remains useful as a future third cross-check, but is no longer required to
close the current Bianchi IX GR_REDUCTION obligation.
