# Repair an existing report summary without rerunning the analysis.
# Open this file in RStudio, edit this ONE path, then run the file.
# On Windows use forward slashes, e.g. "D:/FibrTraQ-results/regional-study".
report_dir <- path.expand("~/FibrTraQ-results/regional-study")

# Requires the FibrTraQ version containing the portable-report-link fix.
# This creates a new summary; the original HTML and scientific files stay intact.
portable_index <- FibrTraQ:::ftq_rebuild_report_index(report_dir)
FibrTraQ:::ftq_qc_open_output(portable_index)
message("Open this summary: ", portable_index)
