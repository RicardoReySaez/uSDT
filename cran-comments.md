# CRAN comments for uSDT 0.1.0

This is a new submission.

## Test environments

* win-builder, R 4.6.1, Windows Server 2022 x64
* win-builder, R-devel 2026-09-09 r90510, Windows Server 2022 x64
* macOS builder, R 4.6.1 Patched 2026-07-27 r90311, macOS 26.6,
  aarch64-apple-darwin23
* Local: Windows 11 x64, R 4.6.0, `R CMD check --as-cran`
* GitHub Actions, defined in `.github/workflows/R-CMD-check.yaml`:
  * macOS-latest, R release
  * windows-latest, R release
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note

The macOS builder reports `Status: OK`. Both win-builder runs report the single
note below, which the macOS builder does not produce because it does not run
the CRAN incoming feasibility check.

### NOTE: checking CRAN incoming feasibility

```
Maintainer: 'Ricardo Rey-Sáez <ricardoreysaez95@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Meyen (26:58)
  SDT (23:57)
  al (26:67)
  et (26:64)
```

This is the first release of uSDT, so the new submission part is expected.

The four flagged words are spelled correctly:

* Meyen is the surname of the first author of the work cited in the
  Description, Meyen et al. (2022) <doi:10.1037/xge0001065>.
* et and al belong to that same citation.
* SDT is the standard abbreviation of signal detection theory. The
  Description writes the term out in full before abbreviating it.

## Notes on the examples

Examples that fit a hierarchical model with lme4 are wrapped in `\donttest{}`
because a single fit takes longer than the five seconds CRAN allows per
example. They all run under `R CMD check --run-donttest`, in about 35 seconds
in total.

One example in `?usdt_boot` is wrapped in `\dontrun{}`. It refits the model 500
times by parametric bootstrap, which takes several minutes on any platform.

## Downstream dependencies

None. This is a new submission, so there are no reverse dependencies to check.
