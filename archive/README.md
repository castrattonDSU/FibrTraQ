# FibrTraQ research archive

The installed package provides `fibrtraq_single()` for one subject and
`fibrtraq_study()` for a study. Both generate the unified report. Ordinary
analyses need only those commands. See the [current workflow guide](../README.md).

This directory preserves the research and validation tools from 0.3.5,
including the in-memory interface, manifest and benchmark scripts, statistical
summaries, historical report tools, and their help and tests.

## Run an earlier research workflow

Install the package from the same checkout. From its repository root, load the
research tools into a separate R environment and call the functions through it.

```r
research <- source("archive/load.R")$value
example <- research$ftq_run_example()
stopifnot(example$validation$passed)
```

The source includes the [previous workflow guide](docs/advanced-workflows-0.3.5.md),
[function help](man/), [research scripts](scripts/) and [validation records](validation/).
Use `research$function_name` for a former package command in an older script.
The archived `research$ftq_run_regional_study()` retains the 0.3.5 default
regional union for local classification. The public `fibrtraq_study()` uses
each focal mask by default.
To reproduce an entire historical run, install the exact revision recorded
with that run and reuse its recorded parameters.

From the repository root, run `Rscript --vanilla archive/tests/run.R` to check
the archived functions against the installed package. The two archived native
audits compile on first use and require an R-compatible C++ compiler. The
optional dTOR plotting script uses `ggplot2` and can use `ggrepel` for labels.
