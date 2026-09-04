# Rebuild the Vadillo et al. (2024) example data.
#
# The package distributes the two original, unmodified data frames as separate
# objects. Keeping them separate mirrors the experiment and lets users pass
# them directly to usdt_data_tasks().

urls <- c(
  awareness = "https://osf.io/bpxe6/download",
  cuing = "https://osf.io/6b2fe/download"
)
expected_md5 <- c(
  awareness = "f62b49d5b8d043a4e998538e307c289c",
  cuing = "f549cf261cf80570f1bc3df1a6f8d952"
)

paths <- vapply(names(urls), function(name) {
  path <- tempfile(fileext = ".csv")
  utils::download.file(urls[[name]], path, mode = "wb", quiet = TRUE)
  path
}, character(1))

actual_md5 <- unname(tools::md5sum(paths))
if (!identical(actual_md5, unname(expected_md5))) {
  stop("The downloaded Vadillo data do not match the archived source files.")
}

vadillo_awareness <- utils::read.csv(paths[["awareness"]])
vadillo_cuing <- utils::read.csv(paths[["cuing"]])

dir.create("data", showWarnings = FALSE)
save(vadillo_awareness, file = "data/vadillo_awareness.rda", compress = "xz")
save(vadillo_cuing, file = "data/vadillo_cuing.rda", compress = "xz")

unlink(paths)
