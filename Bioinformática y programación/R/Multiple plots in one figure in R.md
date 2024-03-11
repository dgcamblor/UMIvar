---
aliases:
  - Combine plots
---

## Base R graphics

Base R provides several functions for creating multiple plots in a single figure. These include:

- `par(mfrow = c(rows, cols))`: This function sets up a grid of rows and columns for creating multiple plots. You can then use standard plotting functions like `plot`, `hist`, or `boxplot` to add plots to this grid. `mfrow` and `mfcol` parameters can be used within the `par()` function to specify the number of rows and columns for arranging plots.
- `layout()`: The `layout` function allows you to specify a custom arrangement of plots in a grid. You define the layout as a matrix, and then you can use the standard plotting functions to fill the individual cells of the grid.

## ggplot2

The options within `ggplot2` to create multiple plots in one figure are the `facet_wrap()` or `facet_grid()` functions, which allow you to split a plot into multiple panels based on one or more variables.

- `facet_wrap()` is used to create a one-dimensional grid of panels based on a single categorical variable.
- `facet_grid()` is used to create a two-dimensional grid of panels based on two categorical variables.

```r
library(ggplot2) 
ggplot(data, aes(x, y)) 
	+ geom_point() 
	+ facet_wrap(~category, ncol = 2)
```

## gridExtra

The `gridExtra` package provides functions like `grid.arrange()` that allow you to arrange multiple plots created with base R graphics or `ggplot2` into a single figure.

Example:

```r
library(gridExtra) 
p1 <- ggplot(data1, aes(x, y)) + geom_point() 
p2 <- ggplot(data2, aes(x, y)) + geom_point() 
grid.arrange(p1, p2, ncol = 2)
```

## cowplot

The `cowplot` package is designed for combining and customizing complex ggplot2 plots. It provides functions like `plot_grid()` for arranging and customizing multiple ggplot2 plots in a single figure.

Example:

```r
library(cowplot)
p1 <- ggplot(data1, aes(x, y)) + geom_point()
p2 <- ggplot(data2, aes(x, y)) + geom_point()
plot_grid(p1, p2, labels = "AUTO")
```