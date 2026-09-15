# Research and publication rules

This repository records exploratory, model-internal numerical experiments. Do not describe a numerical benchmark as an observation, a new physical law, a singularity resolution, or verified quantum gravity.

For each change, state the model, units, initial data, retained/omitted variables, numerical settings, reproducible command, and limitations. Separate implemented results from planned work. Preserve negative controls and failures. Report pair separation separately from forecast error. Do not compare different forecast horizons as equal-horizon accuracy gains.

Before publishing scientific changes, run the relevant tests and record exactly what was rerun. The baseline command is `python reproduce.py`. Optional neural checkpoints are restored from the original ZIP; `python verify_predictors.py` is not a clean-clone prerequisite.

Keep the existing LICENSE. Do not force-push or overwrite unrelated work. Commit source, protocol, results, and interpretation together. Avoid secrets, environment files, and personal conversation transcripts. Update README and the roadmap when the scope changes. No background publication or scheduled execution is configured by this file.
