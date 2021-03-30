### Code for re-generating POMS R package files whenever there has been additional functions added or 
### changes to man pages.

library(devtools)
library(roxygen2)

setwd("~/github_repos/POMS/")
devtools::document()
devtools::load_all(path = "~/github_repos/POMS/")
