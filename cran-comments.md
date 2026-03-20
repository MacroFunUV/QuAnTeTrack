## R CMD check results (Resubmision)
References to methods implemented in the package have been added to the `Description` field in the required CRAN format, using `<doi:...>`, `<ISBN:...>`, and `<https:...>` for auto-linking (as per: https://contributor.r-project.org/cran-cookbook/description_issues.html#references).

Al `print()`, `cat()`, and `writeLines()` statements have been replaced with `message()`, `warning()`, or wrapped in `if (verbose)` where appropriate, following best practices described at: https://contributor.r-project.org/cran-cookbook/code_issues.html#using-printcat.

There is one NOTE about possibly misspelled words ("Tetrapod", "Trackways", "palaeo", "tetrapod", "trackways", author names like "Batschelet", "Benhamou", "Rohlf", etc.); these are domain-specific and intentionally used.

The package now passes `R CMD check --as-cran` with no ERRORs or WARNINGs on my local system.

Thanks for your time and consideration.
