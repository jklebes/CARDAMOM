# CARDAMOM Codebase Metrics

Scope: R post-processing code in `R_functions/` and the Fortran core in
`LIBRARY/CARDAMOM_F/`.

## R code (`R_functions/`, `*.r` / `*.R`)

| # | Metric | Value |
|---|--------|-------|
| 1 | Total lines | 29,224 |
| 2 | Total files | 67 |
| 3 | Total functions | 109 named (111 incl. anonymous) |

## Fortran code (`LIBRARY/CARDAMOM_F/`, `*.f90`)

| # | Metric | Value |
|---|--------|-------|
| 4 | Total lines | 274,203 |
| 5 | Total files | 224 |
| 6 | Functions + subroutines | 2,400 definitions (2,421 incl. interface declarations) |

## Counting methodology

- **R functions** — counted via the named-assignment-of-function pattern
  (`name <- function(`, `<<-`, or `=`): 109. Including the 2 anonymous functions
  passed inline (e.g. to `apply`), the total `function(` occurrences are 111.
- **Fortran procedures** — the headline 2,400 counts true procedure *definitions*
  (1,765 subroutines + 635 functions), excluding `interface`-block declarations.
  The raw `end subroutine`/`end function` tally is 2,421 (1,772 + 649); the
  21-procedure difference is the interface-block prototype declarations (which
  legitimately end with `end subroutine`/`end function` but are not
  implementations).
