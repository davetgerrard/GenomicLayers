
# runs all the scripts in the tests/ sub-folder. Reports on failed tests.


library(testthat) 

source('scripts/predictFromSequence.functions.R')

test_results <- test_dir("tests", reporter="summary")



