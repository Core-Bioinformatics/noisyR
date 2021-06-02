.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to noisyr!")
  packageStartupMessage("Have a look at our website to get started:")
  packageStartupMessage("  https://core-bioinformatics.github.io/noisyR/")
  packageStartupMessage("To cite noisyr in publications please use this paper:")
  packageStartupMessage(paste0("  ", citation("noisyr")[[1]]$url))
}
