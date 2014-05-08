hab
===

This is the development area for the package `hab`, a set of functions
related to habitat selection and movement analyses. Also includes a few
patches for functions from the adehabitatXY package family.

See [the home page for the released versions](http://ase-research.org/basille/hab/).


## Installation of the development version

You need to use the package
[`devtools`](http://cran.r-project.org/web/packages/devtools/index.html)
from Hadley Wickham:

	library(devtools)
	install_github("basille/hab")


## Installation of the released versions

Binaries for Windows are (sometimes) maintained for the
latest R release, and can be installed with the command:

    install.packages("hab", repos = "http://ase-research.org/R/")

If it's not the case, and with other operating systems (such as GNU/Linux or OS X), you can compile the source package with the command:

	install.packages("hab", repos = "http://ase-research.org/R/", type = "source")

The upgrade process (e.g. using `update.packages()`) is not operative on
packages from this repository. Simply reinstall the package with one of the
commands above to upgrade it when necessary.


## Getting started

For a list of documented functions, use `library(help = "hab")` or see the
Reference manual.
