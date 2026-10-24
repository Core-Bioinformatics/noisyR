## Load the noisyR functions from the package sources in this checkout, without
## installing the package: every R/*.R file is sourced into one environment with
## the `noisyr::` prefixes stripped, and `%>%` taken from magrittr.
## Run the scripts from the repository root, e.g.
##   R_LIBS=<lib with philentropy, preprocessCore, ...> Rscript analysis/grn_noise/run_benchmark.R
load_noisyr <- function(pkg.dir = "."){
  env <- new.env()
  assign("%>%", magrittr::`%>%`, envir = env)
  files <- list.files(file.path(pkg.dir, "R"), pattern = "\\.R$", full.names = TRUE)
  files <- files[!(basename(files) %in% c("zzz.R", "run_noisyr_app.R"))]
  for(f in files){
    code <- gsub("noisyr::", "", readLines(f))
    eval(parse(text = code), envir = env)
  }
  env
}
