/-
Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0

Finite rational arithmetic certificates, not a formalization of real integration
or of the Einstein-Vlasov evolution. Coefficient k represents x^(2*k).
-/
import Std.Internal.Rat

namespace AstraBlackhole
open Std.Internal

/-- Finite functional: sum_k c[k] / (2*(shift+k)+1). -/
def evenMoment : List Rat → Nat → Rat
  | [], _ => 0
  | c :: cs, shift =>
      c * mkRat 1 (2 * shift + 1) + evenMoment cs (shift + 1)

-- Even-power coefficients of P4, P6, P8, in ascending order.
def p4 : List Rat := [3/8, -30/8, 35/8]
def p6 : List Rat := [-5/16, 105/16, -315/16, 231/16]
def p8 : List Rat := [35/128, -1260/128, 6930/128, -12012/128, 6435/128]

theorem p4_low_moments :
    evenMoment p4 0 = 0 ∧ evenMoment p4 1 = 0 := by decide

theorem p4_first_visible : evenMoment p4 2 = 8/315 := by decide

theorem p6_low_moments :
    evenMoment p6 0 = 0 ∧ evenMoment p6 1 = 0 ∧
    evenMoment p6 2 = 0 := by decide

theorem p6_first_visible : evenMoment p6 3 = 16/3003 := by decide

theorem p8_low_moments :
    evenMoment p8 0 = 0 ∧ evenMoment p8 1 = 0 ∧
    evenMoment p8 2 = 0 ∧ evenMoment p8 3 = 0 := by decide

theorem p8_first_visible : evenMoment p8 4 = 128/109395 := by decide

end AstraBlackhole
