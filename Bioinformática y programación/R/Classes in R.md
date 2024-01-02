There are three main systems for working with classes in R: S3, S4 and S6 (Reference Class).

## S3

S3 classes are a simple, straightforward system of creating classes in R. It's often more convenient than S4 for quick tasks.

```R
# Create an S3 object with a class attribute
genetic_variant <- list(
  CHR = "chr1",
  POSITION = 1000,
  REF = "A",
  ALT = "G",
  VAF = 0.05
)

# Assign a class attribute to the object
class(genetic_variant) <- "GeneticVariant"

# Define a generic function
get_CHR <- function(variant) {
  if (class(variant) == "GeneticVariant") {
    return(variant$CHR)
  } else {
    stop("Invalid object class.")
  }
}

# Use the generic function
chr_value <- get_CHR(genetic_variant)
cat("CHR:", chr_value, "\n")

```

In S3 systems, methods don’t belong to the class, but instead belong to generic functions.
## S4

S4 systems are a more formal way to define classes. The structure of S4 classes is more akin to that of other languages like C++ or Java.

```r
setClass(
  "GeneticVariant",
  representation(
    CHR = "character",
    POSITION = "numeric",
    REF = "character",
    ALT = "character",
    VAF = "numeric",
    TYPE = "character"
  )
)

get_variant_type <- function(ref, alt) {
  if (nchar(ref) == 1 & nchar(alt) == 1) {
    return("SNV")
  } else if (nchar(ref) > 1 & nchar(alt) > 1) {
    return("MNV")
  } else if (nchar(ref) == 1 & nchar(alt) > 1) {
    return("DEL")
  } else if (nchar(ref) > 1 & nchar(alt) == 1) {
    return("INS")
  } else {
    return("COMPLEX")
  }
}

setMethod(
  "initialize",
  "GeneticVariant",
  function(.Object, CHR, POSITION, REF, ALT, VAF) {
    .Object@CHR <- CHR
    .Object@POSITION <- POSITION
    .Object@REF <- REF
    .Object@ALT <- ALT
    .Object@VAF <- VAF
    .Object@TYPE <- get_variant_type(REF, ALT)
    return(.Object)
  }
)

setMethod(
  "show",
  "GeneticVariant",
  function(object) {
    cat(paste0(object@CHR, ":", object@POSITION, object@REF, ">", object@ALT, "(", object@VAF, ")"))
  }
)

# Getters
setGeneric("get_variant", function(object) {standardGeneric("get_variant")})

setMethod(
  "get_variant",
  "GeneticVariant",
  function(object) {
    return(paste0(object@CHR, ":", object@POSITION, object@REF, ">", object@ALT))
  }
)
```