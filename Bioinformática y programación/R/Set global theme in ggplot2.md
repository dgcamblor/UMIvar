---
tags:
  - how-to
---
```r
# Set global theme options
my_theme <- theme_minimal() +
  theme(
    text = element_text(color = "blue"),  # Set text color to blue
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  )
```