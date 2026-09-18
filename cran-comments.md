# CRAN comments for uSDT 0.1.0

## Resubmission

This is a resubmission of a new package. In this version I have:

* Replaced `\dontrun{}` with `\donttest{}` in the example of `usdt_boot()`.
* Unwrapped every example that runs in less than five seconds. Only the two
  blocks described in "Notes on the examples" remain in `\donttest{}`.
* Removed the code in `usdt_boot()` that created `.Random.seed` in the global
  environment when it did not exist. The function now changes the random seed
  only through `set.seed()`, when the user supplies `seed`, which defaults to
  `NULL`.
* Changed the order of the authors in `Authors@R`.

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

The examples of each help page run in less than five seconds, except two
blocks that stay in `\donttest{}`:

* The second part of `?usdt_hypotheses` fits a trial-level model with
  `lme4::glmer()`, which takes about 25 seconds.
* `?usdt_boot` refits the model 500 times by parametric bootstrap, which takes
  about three minutes. 500 is the fewest replicates `usdt_boot()` accepts.

Both blocks run under `R CMD check --run-donttest`.

## Downstream dependencies

None. This is a new submission, so there are no reverse dependencies to check.
