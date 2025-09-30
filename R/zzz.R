# Ensure g.anat is the canonical object when the namespace loads (CI-safe)
.onLoad <- function(libname, pkgname) {
  assign("gganat", .build_gganat(), envir = asNamespace(pkgname))
}

