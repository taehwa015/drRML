library(knitr)
# remove file path vignettes
replace = readLines("vignettes/drRML.Rmd")
replace = gsub("<img src=\"vignettes/", "<img src=\"", replace)
fileConn = file("vignettes/drRML.Rmd")
writeLines(replace, fileConn)
close(fileConn)
