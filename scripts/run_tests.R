
# runs all the scripts in the tests/ sub-folder. Reports on failed tests.


library(testthat) 

source('scripts/predictFromSequence.functions.R')

test_results <- test_dir("tests/testthat", reporter="summary")


write.table(test_results, file="testResults.tsv", quote=F, row.names=F, sep="\t")
