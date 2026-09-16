# Research and publication rules

This repository records exploratory, model-internal numerical experiments. Do not describe a numerical benchmark as an observation, a new physical law, a singularity resolution, or verified quantum gravity.

For each change, state the model, units, initial data, retained/omitted variables, numerical settings, reproducible command, and limitations. Separate implemented results from planned work. Preserve negative controls and failures. Report pair separation separately from forecast error. Do not compare different forecast horizons as equal-horizon accuracy gains.

Before publishing scientific changes, run the relevant tests and record exactly what was rerun. The baseline command is `python reproduce.py`. Optional neural checkpoints are restored from the original ZIP; `python verify_predictors.py` is not a clean-clone prerequisite.

Preserve the Apache-2.0 LICENSE and NOTICE (Copyright 2026 HeliCorgi), unless the owner explicitly authorizes another change. The change from MIT was explicitly authorized on 2026-09-16. Retain historical provenance and any third-party notices. Do not force-push or overwrite unrelated work. Commit source, protocol, results, and interpretation together. Avoid secrets, environment files, and personal conversation transcripts. Update README and the roadmap when the scope changes.

For Lean changes, run `cd lean && lake build && lake env lean -DwarningAsError=true Audit.lean`, check the axiom report, and record whether compilation actually completed. Do not use placeholders or add axioms to make a claim pass. Keep finite arithmetic certificates separate from real-integral identities, evolution theorems, and physical validity. A configured workflow or a queued run is not a passed proof.

CI runs only on push, pull request, and manual dispatch. No scheduled research or automatic repository-writing execution is configured by these files.
