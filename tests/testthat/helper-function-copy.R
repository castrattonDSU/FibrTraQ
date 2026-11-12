# Intercept selected export lookups only within a copied function. The original
# function and every package namespace retain their bindings and environments.
ftq_test_function_copy <- function(fun, exports = list(), bindings = list()) {
  lookup_export <- function(package, name) {
    package <- as.character(substitute(package))
    name <- as.character(substitute(name))
    key <- paste0(package, "::", name)
    if (key %in% names(exports)) return(exports[[key]])
    getExportedValue(package, name)
  }
  environment(fun) <- list2env(
    c(bindings, list(`::` = lookup_export)), parent = environment(fun)
  )
  fun
}
