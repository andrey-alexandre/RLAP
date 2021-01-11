library(devtools)
library(testthat)
library(tidyverse)
library(fs)

use_git()

a <- factor(c("character", "hits", "your", "eyeballs"))
b <- factor(c("but", "integer", "where it", "counts"))

use_r("fbind")
load_all()

fbind(a, b)

check()

use_mit_license("Andrey Alexandre")

document()

check()
install()

use_testthat()

use_test("fbind")

test()

use_package("forcats")

use_r("fcount")

document()

use_readme_rmd()
