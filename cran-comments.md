# CRAN comments for uSDT 0.1.0

This is a new submission.

## Test environments

* Local: Windows 11 x64, R 4.6.0, `R CMD check --as-cran`
* GitHub Actions, defined in `.github/workflows/R-CMD-check.yaml`:
  * macOS-latest, R release
  * windows-latest, R release
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1

<!-- Before submitting:
     1. Confirm the GitHub Actions matrix above is green for the commit
        being submitted.
     2. Run devtools::check_win_devel(), devtools::check_win_release() and
        devtools::check_mac_release(), and add their results here.
     3. Delete this comment block. -->

## R CMD check results

0 errors | 0 warnings | 2 notes

### NOTE: checking CRAN incoming feasibility

```
Maintainer: 'Ricardo Rey-Sáez <ricardoreysaez95@gmail.com>'

New submission
```

This is the first release of uSDT, so the note is expected.

### NOTE: checking HTML version of manual

```
Skipping checking math rendering: package 'V8' unavailable
```

The local check machine does not have V8 installed. The note does not appear
where V8 is available.

## Notes on the examples

Examples that fit a hierarchical model with lme4 are wrapped in `\donttest{}`
because a single fit takes longer than the five seconds CRAN allows per
example. They all run under `R CMD check --run-donttest`, in about 35 seconds
in total.

One example in `?usdt_boot` is wrapped in `\dontrun{}`. It refits the model 500
times by parametric bootstrap, which takes several minutes on any platform.

## Downstream dependencies

None. This is a new submission, so there are no reverse dependencies to check.
