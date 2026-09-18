/-
Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0

Semantic audit lemmas only. These theorems do not choose a physical quantum
constraint, clock, inner product, operator domain, factor ordering, or regulator.
-/
import Init

namespace AstraBlackhole

/-- If a first-order branch equation, an abstract product-rule correction, and
an abstract square relation hold, the corresponding second-order equation
contains exactly that residual term. The operation `add` is left abstract. -/
theorem time_dependent_factorization_residual
    {State : Type}
    (Ps G A R : State → State)
    (add : State → State → State)
    (psi : State)
    (hbranch : Ps psi = G psi)
    (hproduct : Ps (G psi) = add (G (G psi)) (R psi))
    (hsquare : G (G psi) = A psi) :
    Ps (Ps psi) = add (A psi) (R psi) := by
  calc
    Ps (Ps psi) = Ps (G psi) := congrArg Ps hbranch
    _ = add (G (G psi)) (R psi) := hproduct
    _ = add (A psi) (R psi) := congrArg (fun x => add x (R psi)) hsquare

/-- Exact branch recombination remains exact after applying the same operator
U to the recombined state. -/
theorem common_operator_preserves_recombination
    {State BranchData : Type}
    (split : State → BranchData)
    (recombine : BranchData → State)
    (U : State → State)
    (psi : State)
    (hrecombine : recombine (split psi) = psi) :
    U (recombine (split psi)) = U psi :=
  congrArg U hrecombine

/-- A zero-preserving operator K that commutes with the constraint C on psi
maps a constraint-kernel state to another constraint-kernel state. -/
theorem commuting_operator_preserves_constraint_kernel
    {State : Type}
    (C K : State → State)
    (zero psi : State)
    (hcommute : C (K psi) = K (C psi))
    (hconstraint : C psi = zero)
    (hzero : K zero = zero) :
    C (K psi) = zero :=
  hcommute.trans ((congrArg K hconstraint).trans hzero)

end AstraBlackhole
