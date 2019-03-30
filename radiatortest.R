
library(dplyr)
# please re-install, i've just made a commit
devtools::install_github("thierrygosselin/radiator")

BiocManager::install("SeqArray")
BiocManager::install("zlibbioc")
#BiocManager::install("gdfsmt")

library("devtools")
install_github("zhengxwen/gdsfmt")

BiocManager::install("SeqVarTools")	

library(radiator)



test.gl <- radiator::read_dart(data = "d:/temp/dart_count_dartr_test.tsv", strata = "d:/temp/dart_strata_dartr_test.tsv") %>%
  radiator::write_genlight(data = ., dartr = TRUE)

test.gl@other$loc.metrics <- data.frame(test.gl@other$loc.metrics)
#checks
names(test.gl@other$loc.metrics)
# I’ll let you test more…
