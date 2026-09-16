/-
Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
-/
import Init

namespace AstraBlackhole

/-- Equal observations yield equal deterministic predictions. -/
theorem equal_observations_equal_predictions
    {State Observation Prediction : Type}
    (observe : State → Observation) (forecast : Observation → Prediction)
    (a b : State) (h : observe a = observe b) :
    forecast (observe a) = forecast (observe b) :=
  congrArg forecast h

/-- A deterministic predictor cannot exactly recover two different futures
from identical observations. No existence of such physical states is asserted. -/
theorem no_exact_forecast_for_distinct_futures
    {State Observation Prediction : Type}
    (observe : State → Observation) (forecast : Observation → Prediction)
    (future : State → Prediction) (a b : State)
    (hobs : observe a = observe b) (hfuture : future a ≠ future b) :
    ¬ (forecast (observe a) = future a ∧ forecast (observe b) = future b) := by
  intro h
  apply hfuture
  exact h.1.symm.trans ((congrArg forecast hobs).trans h.2)

end AstraBlackhole
