script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_file),
    winslash = "/",
    mustWork = TRUE
  )
  ROOT <- dirname(dirname(script_path))
} else {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

project_library <- file.path(ROOT, ".Rlib")
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_library, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = project_library)
}

bioconductor_packages <- c("edgeR", "limma", "Rsubread", "tximport")
missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_bioconductor)) {
  BiocManager::install(
    missing_bioconductor,
    lib = project_library,
    ask = FALSE,
    update = FALSE
  )
}

cran_packages <- c("readxl", "digest", "jsonlite")
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = project_library)
}

cat(
  "R packages ready in",
  normalizePath(project_library, winslash = "/"),
  "\n"
)
