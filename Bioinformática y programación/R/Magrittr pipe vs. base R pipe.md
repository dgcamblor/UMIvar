
By default, the `|>` pipe passes the object on its left-hand side as the first argument, similar to `%>%`. Some of the key differences are:

- `%>%` allows changing the placement with a `.`, while `|>` uses `_` as a placeholder with the additional restriction that the argument has to be named.
- `|>` placeholder is simpler, lacking some features of `%>%`, such as passing to multiple arguments or having special behavior inside another function.
- `%>%` supports using `.` on the left-hand side of operators like `$`, `[[`, `[`, enabling operations like extracting a single column from a data frame with `mtcars %>% .$cyl`.
- `%>%` allows dropping parentheses when calling a function with no other arguments, whereas `|>` always requires parentheses.
- Because the native pipe was introduced in R 4.1.0, using `|>` in function reference examples or vignettes may not work on older R versions.

## References

[Differences between the base R and magrittr pipes](https://www.tidyverse.org/blog/2023/04/base-vs-magrittr-pipe/)