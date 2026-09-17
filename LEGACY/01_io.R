
# 01_io.R — reading data

read_bed_like <- function(path, sep = ";") {
  read.table(
    file = path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
