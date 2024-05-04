# Enterobase webpage scraper


## Description:
This R script allows
   - listing all of the available schemes on [Enterobase](https://enterobase.warwick.ac.uk/schemes/)
   - listing all of the available organisms on Enterobase
   - finding the date and time of the last change on scheme of interest
   - downloading schemes from Enterobase


## Requirements
- [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
- [rvest](https://cran.r-project.org/web/packages/rvest/index.html)
- [knitr](https://cran.r-project.org/web/packages/knitr/index.html)
- [optparse](https://cran.r-project.org/web/packages/optparse/index.html)

In case you would like to use this script you can easily install all the required packages by running the code below in your R session: 

```R
# Listing required packages
required_packages <- c("tidyverse", "rvest", "knitr", "optparse")

# Check if required packages are installed
missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])

# Install missing packages
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}
```


## Usage:

To see help

```bash
Rscript --vanilla enterobase_scheme_scraper.R --help
```

To list available organisms on Enterobase

```bash
Rscript --vanilla enterobase_scheme_scraper.R -f list_organisms
```

To list available schemes for given organism

```bash
Rscript --vanilla enterobase_scheme_scraper.R -f list_organism_schemes -o Salmonella
```

To download scheme

```bash
Rscript --vanilla enterobase_scheme_scraper.R -f download_scheme -o Salmonella -s Achtman7GeneMLST
```
