################################################################################
# Enterobase webpage scraper
#
# Author: Vladimir Bajić
# Date: 2024-05-04
#
# Description:
# This script allows
#   - listing all of the available schemes on Enterobase (https://enterobase.warwick.ac.uk/schemes/)
#   - listing all of the available organisms on Enterobase
#   - finding the date and time of the last change on scheme of interest
#   - downloading schemes from Enterobase
#
#
# Usage:
#
# To see help
#   Rscript --vanilla enterobase_scheme_scraper.R --help
#
# To list available organisms on Enterobase
#   Rscript --vanilla enterobase_scheme_scraper.R -f list_organisms
#
# To list available schemes for given organism
#   Rscript --vanilla enterobase_scheme_scraper.R -f list_organism_schemes -o Salmonella
#
# To download scheme
#   Rscript --vanilla enterobase_scheme_scraper.R -f download_scheme -o Salmonella -s Achtman7GeneMLST
#
################################################################################

# Libraries --------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(rvest))
library(knitr)
library(optparse)

# Define base URL --------------------------------------------------------------
url <- "https://enterobase.warwick.ac.uk/schemes/"

# Functions --------------------------------------------------------------------

## Function to get scheme index table from url
scheme_index_table_from_url <- function(url) {
    ### Save a nodes from html
    html_nodes_a <-
        url %>%
        read_html() %>%
        html_nodes("a")

    ### Remove nodes that contain href="../" or href="%2A.%2A"
    html_nodes_a_clean <-
        html_nodes_a[
            !grepl("\\.\\.", html_attr(html_nodes_a, "href")) &
                !grepl("%2A.%2A", html_attr(html_nodes_a, "href"))
        ]

    ### Extract table with info for each scheme
    scheme_tbl <-
        html_nodes_a_clean %>%
        html_attr("href") %>%
        as_tibble_col(column_name = "URL_path") %>%
        separate(col = URL_path, into = c("Organism", "Scheme"), sep = "\\.", remove = FALSE) %>%
        mutate(Scheme_ID = str_remove(URL_path, "/")) %>%
        mutate(Scheme = str_remove(Scheme, "/")) %>%
        mutate(Full_Path = paste0(url, URL_path))

    ### Extract information with dates and times of last change
    dt_tbl <-
        html_nodes_a_clean %>%
        html_nodes(xpath = "./following-sibling::text()") %>%
        html_text2() %>%
        str_remove(., "\r") %>%
        as_tibble_col(column_name = "Date_Time_Size") %>%
        separate(col = Date_Time_Size, into = c("Date", "Time", "Size"), sep = " ") %>%
        mutate(LastUpdated = paste(Date, Time, sep = "_")) %>%
        select(-Size)

    ### Bind columns containing dates and table info into one table
    scheme_index_table <- bind_cols(scheme_tbl, dt_tbl)

    return(scheme_index_table)
}

## Function to print formatted list of all schemes from Enterobase and their info
list_schemes <- function() {
    scheme_index_table %>%
        select(Organism, Scheme, Date, Time) %>%
        kable() %>%
        print()
}

## Function to list all available organisms on pubMLST
list_organisms <- function() {
    scheme_index_table %>%
        select(Organism) %>%
        distinct() %>%
        kable() %>%
        print()
}

## Function to list available schemes for given organism_id unformatted
list_organism_schemes <- function(organism_id) {
    scheme_index_table %>%
        filter(Organism == organism_id) %>%
        select(Organism, Scheme, Date, Time) %>%
        distinct() %>%
        kable() %>%
        print()
}

## Function to download fasta files of alleles included in schema
download_scheme <- function(oid, sid) {
    ### Select only relevant row
    scheme_index_table_sub <-
        scheme_index_table %>%
        filter(Organism == oid & Scheme == sid)

    ### Save html
    tmp_html <-
        scheme_index_table %>%
        filter(Organism == oid & Scheme == sid) %>%
        pull(Full_Path)

    ### Save a nodes from html
    html_nodes_a <-
        tmp_html %>%
        read_html() %>%
        html_nodes("a")

    ### Remove nodes that contain href="../" or href="%2A.%2A"
    html_nodes_a_clean <-
        html_nodes_a[
            !grepl("\\.\\.", html_attr(html_nodes_a, "href")) &
                !grepl("%2A.%2A", html_attr(html_nodes_a, "href"))
        ]

    ### Extract table to download
    file_tbl <-
        html_nodes_a_clean %>%
        html_attr("href") %>%
        as_tibble_col(column_name = "File") %>%
        mutate(Download = paste0(tmp_html, File))

    ### Extract information with dates and times of last change
    dts_tbl <-
        html_nodes_a_clean %>%
        html_nodes(xpath = "./following-sibling::text()") %>%
        html_text2() %>%
        str_remove(., "\r") %>%
        as_tibble_col(column_name = "Date_Time_Size") %>%
        separate(col = Date_Time_Size, into = c("Date", "Time", "Size"), sep = " ")

    ### Bind columns containing dates and table info into one table
    file_info_tbl <- bind_cols(file_tbl, dts_tbl)

    ### Make scheme timestamp for the name of the scheme to be downloaded
    scheme_timestamp <- paste0("schemeID_", scheme_index_table_sub$Scheme_ID, "_LastUpdated_", scheme_index_table_sub$LastUpdated)

    ### Create path to where the scheme will be downloaded
    destfile_path <- paste0(scheme_timestamp)

    ### Check if dir already exists and if yes do not download
    if (dir.exists(destfile_path)) {
        stop("WARNING: ", destfile_path, " already exists. \nDownloading aborted to prevent overwriting.\n")
    }

    ### Create dir where all files will be downloaded
    dir.create(destfile_path)
    dir.create(paste0(destfile_path, "/loci_fastas/"))


    ### Open a connection to a log file in append mode
    log_file <- file(paste0(destfile_path, "/download_error.log"), open = "a")

    ### Initialize counters
    warning_count <- 0
    error_count <- 0

    ### Loop to download all alleles
    for (i in seq_along(file_info_tbl$Download)) {
        cat("Downloading File", i, "/", length(file_info_tbl$Download), "\n")

        tryCatch(
            {
                ##### Check if the file is fasta.gz and if so save it in special sub directory
                if (grepl(".fa.*.gz$", file_info_tbl$Download[i])) {
                    download.file(
                        file_info_tbl$Download[i],
                        paste0(destfile_path, "/loci_fastas/", file_info_tbl$File[i], sep = ""),
                        mode = "auto", quiet = TRUE
                    )
                } else {
                    download.file(
                        file_info_tbl$Download[i],
                        paste0(destfile_path, "/", file_info_tbl$File[i], sep = ""),
                        mode = "auto", quiet = TRUE
                    )
                }
            },
            error = function(e) {
                # Increment error count
                error_count <<- error_count + 1
                # Write error message to log file
                cat("Error downloading:", basename(file_info_tbl$Download[i]), "-", conditionMessage(e), "\n")
                writeLines(paste("Error downloading:", basename(file_info_tbl$Download[i]), "-", conditionMessage(e)), log_file)
            },
            warning = function(w) {
                # Increment warning count
                warning_count <<- warning_count + 1
                # Write warning message to log file
                cat("Warning downloading:", basename(file_info_tbl$Download[i]), "-", conditionMessage(w), "\n")
                writeLines(paste("Warning downloading:", basename(file_info_tbl$Download[i]), "-", conditionMessage(w)), log_file)
            }
        )
    }

    # Close the log file connection
    close(log_file)

    # Output counts
    if (warning_count + error_count > 0) {
        cat("\nWARNING: Number of files that were not downloaded:", warning_count + error_count, "\n")
    }
}

#-------------------------------------------------------------------------------

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-f", "--function"),
        type = "character", metavar = "character",
        help = "Function to be performed \n\n\t\tPossible functions are:\n
            \t| list_schemes              |   to list all available schemes on Enterobase
            \t| list_organisms            |   to list all available organisms on Enterobase
            \t| list_organism_schemes     |   to list all available Enterobase schemes for a given organism_id
            \t| download_scheme           |   to download all files included in schema \n"
    ),
    make_option(c("-o", "--organismID"),
        type = "character", metavar = "character",
        help = "Organism ID on which to perform function\n"
    ),
    make_option(c("-s", "--schemeID"),
        type = "character", metavar = "character",
        help = "Scheme ID on which to perform function\n"
    )
)
# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Make scheme index table from URL00 -------------------------------------------
scheme_index_table <- scheme_index_table_from_url(url)

# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$f)) {
    print_help(opt_parser)
    stop("Choose one of the possible functions.")
}

if (opt$f == "list_schemes") {
    list_schemes()
}

if (opt$f == "list_organisms") {
    list_organisms()
}

if (opt$f == "list_organism_schemes") {
    if (is.null(opt$o)) {
        cat("\nPlease provide Organism_ID using argument '-o' for which you want to see the list of available schemes\n\n")
    } else {
        list_organism_schemes(opt$o)
    }
}

if (opt$f == "download_scheme") {
    if (is.null(opt$o) | is.null(opt$s)) {
        cat("\nPlease provide Organism ID and Scheme ID using arguments '-o' and '-s'.\n\n")
    } else {
        download_scheme(opt$o, opt$s)
    }
}
