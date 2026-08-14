# CRAN comments for zzrctsim 0.1.0

## Submission

This is a new submission.

## Test environments

- local: macOS 26.5.1 (aarch64-apple-darwin25.4.0), R 4.6.1
- (to be added before submission: win-builder devel and release,
  R-hub linux/windows/macos)

## R CMD check results

0 errors | 0 warnings | 1 note

The one NOTE is:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ronald G. Thomas <rgthomas@ucsd.edu>'
New submission
```

This is expected for a first submission.

## Notes for the reviewer

The DESCRIPTION references two methods papers, both cited with DOIs in
the standard `<doi:...>` form: Diggle and Kenward (1994)
<doi:10.2307/2986113> for the selection-model dropout mechanism, and
Morris, White and Crowther (2019) <doi:10.1002/sim.8086> for the
simulation-study reporting standard the package follows.

Examples and tests run quickly; the simulation examples use small
replicate counts so that no example exceeds the CRAN time limits.
