# Ensure g.anat is the canonical object when the namespace loads (CI-safe)
.onLoad <- function(libname, pkgname) {
  assign("g.anat", .build_g_anat(), envir = asNamespace(pkgname))
}

