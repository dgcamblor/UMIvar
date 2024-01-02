---
url: 
tags:
  - package
---

To create publication-quality TIFF Venn Diagrams:

```r
venn.diagram(
  x = list(inc_tps, mageri_tps, uvc_tps), # Vectors of elements
  category.names = c("INCLIVA", "MAGERI", "UMI-VarCal"),
  filename = "venn.tiff",
  # Appearance
  fill = inc_palette[1:3],
  col = inc_palette_darker[1:3],
  lwd = 3,
  cex = 1.3,
  cat.dist = 0.05,
  fontface = "bold",
  cat.fontface = "bold"
)
# Other figure arguments: height, width, resolution (DPI)
```

The `venn.diagram()` function operates in two modes:

- If `filename = "filename.tiff"`, the output of the function is a file.
- If `filename = NULL`, the output is an image object, which can be passed to the `grid.draw()` function to produce a plot.