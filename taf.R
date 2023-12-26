#install.packages("icesTAF")

library(icesTAF)
run_dir <- download.analysis("ices-taf-dev/wiki-example-stockassessment")
# view downloaded code
browseURL(run_dir)

run_dir <- download.analysis("ices-taf/2019_san.sa.6", dir = ".")
browseURL(run_dir)

rm(run_dir)
run_dir <- download.analysis("ices-taf/2022_bss.27.4bc7ad-h_assessment")   #public
browseURL(run_dir)

rm(run_dir)
run_dir <- download.analysis("ices-taf/2022_cod.27.47d20_assessment") # not public (as most of the assessments)
browseURL(run_dir)
