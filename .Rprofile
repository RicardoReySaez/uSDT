# .Rprofile
# This file prepares the R sessions that are started in the package directory.
# Author: Ricardo Rey-Sáez
# Last modified: 05-09-2026

# R reads a single profile, and a file here takes precedence over the one in
# the home directory, so without this line a project profile would silently
# replace the user's own settings rather than add to them.
if (file.exists("~/.Rprofile")) source("~/.Rprofile")

# The print methods draw their rules with box-drawing characters and write the
# group difference as a Greek delta. R transliterates both on the way out when
# LC_CTYPE is not UTF-8: every rule becomes a row of hyphens and the delta
# becomes "Delta". A site built from such a session then records the locale of
# whoever ran the build rather than what the package printed, in the reference
# examples as much as in the vignette, because the loss happens when the output
# is written and not when the characters are chosen.
#
# The request is allowed to fail. A machine without that locale keeps the one
# it had and everything still runs; only those characters degrade.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
}
