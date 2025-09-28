pkgname <- "ubynet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ubynet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("build_mass_pmd_network")
### * build_mass_pmd_network

flush(stderr()); flush(stdout())

### Name: build_mass_pmd_network
### Title: Build Mass-PMD Transformation Network
### Aliases: build_mass_pmd_network

### ** Examples

## Not run: 
##D build_mass_pmd_network(
##D     mol_file = "MS_MolInfor1.csv",
##D     trans_file = "Transformation_Database_07-2020.csv",
##D     error_term = 0.00001,
##D     output_dir = "MS_MolInfor2"
##D )
## End(Not run)




cleanEx()
nameEx("classify_Mass")
### * classify_Mass

flush(stderr()); flush(stdout())

### Name: classify_Mass
### Title: Classify Mass into product, resistant, and resistant
### Aliases: classify_Mass

### ** Examples

# Using data frames
before_data <- data.frame(Mass = c(180.16, 46.07))
after_data <- data.frame(Mass = c(46.07, 16.04))
classify_Mass(before_data, after_data)

# Using CSV files
classify_Mass("before.csv", "after.csv", "classified_results.csv")



cleanEx()
nameEx("classify_MolForm")
### * classify_MolForm

flush(stderr()); flush(stdout())

### Name: classify_MolForm
### Title: Classify MolForm into product, resistant, and disappearance
### Aliases: classify_MolForm

### ** Examples

# Using data frames
before_data <- data.frame(MolForm = c("C6H12O6", "C2H5OH"))
after_data <- data.frame(MolForm = c("C2H5OH", "CH4"))
classify_MolForm(before_data, after_data)

# Using CSV files
classify_MolForm("before.csv", "after.csv", "classified_results.csv")



cleanEx()
nameEx("compare_mass")
### * compare_mass

flush(stderr()); flush(stdout())

### Name: compare_mass
### Title: Compare Mass columns between two CSV files
### Aliases: compare_mass

### ** Examples

## Not run: 
##D result <- compare_mass("MS_MolInfor1.csv", "MS_MolInfor2.csv", 
##D                         output_dir = "results", mass_tolerance = 0.01)
## End(Not run)



cleanEx()
nameEx("compare_molforms")
### * compare_molforms

flush(stderr()); flush(stdout())

### Name: compare_molforms
### Title: Compare MolForm columns between two CSV files
### Aliases: compare_molforms

### ** Examples

## Not run: 
##D result <- compare_molforms("MS_MolInfor1.csv", "MS_MolInfor2.csv", output_dir = "results")
## End(Not run)



cleanEx()
nameEx("compare_multiple_datasets")
### * compare_multiple_datasets

flush(stderr()); flush(stdout())

### Name: compare_multiple_datasets
### Title: Compare multiple CSV datasets
### Aliases: compare_multiple_datasets

### ** Examples

## Not run: 
##D files <- c("MS_MolInfor1.csv", "MS_MolInfor2.csv", "MS_MolInfor3.csv")
##D results <- compare_multiple_datasets(files, output_dir = "results", comparison_type = "molform")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
