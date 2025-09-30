# Ensure the canonical gganat object is present when the namespace loads.
.onLoad <- function(libname, pkgname) {
  assign("gganat", .build_gganat(), envir = asNamespace(pkgname))
}
