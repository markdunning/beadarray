
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Welcome to beadarray version ", packageDescription("beadarray", field="Version"))

  packageStartupMessage("beadarray versions >= 2.0.0 are substantial updates from beadarray 1.16.0 and earlier. Please see package vignette for details")

}

