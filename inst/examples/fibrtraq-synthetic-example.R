# Reproduce the bundled FibrTraQ known-answer dataset from a fresh R session.

library(FibrTraQ)

example <- ftq_run_example(
  output_dir = file.path(tempdir(), "fibrtraq-known-answer-output"),
  overwrite = TRUE,
  verbose = FALSE
)

print(example)
print(example$validation$tables)

cat("\nStandard report:    ", example$primary_report, "\n", sep = "")
cat("Exploratory report: ", example$exploratory_report, "\n", sep = "")
cat("Interactive 3-D QC: ", example$interactive_3d_qc, "\n", sep = "")

invisible(example)
