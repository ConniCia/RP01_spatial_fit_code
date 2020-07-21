# Paste with underscore separator ----------------------------------------------
## remove any "_" right before a "."
## reduce multiple "_" in sequence to a single "_"
## remove leading and trailing "_"

paste_ <- function(...) {
  str <- paste(..., sep = "_")
  str <- gsub("_[.]", ".", str)
  str <- gsub("_{2,}", "_", str)
  str <- gsub("^_", "", str)
  gsub("_$", "", str)
}


# Paste with fullstop separator ------------------------------------------------
## reduce multiple "." in sequence to a single "."
## remove leading and trailing "."

paste. <- function(...) {
  str <- paste(..., sep = ".")
  str <- gsub("[.]{2,}", ".", str)
  str <- gsub("^[.]", "", str)
  gsub("[.]$", "", str)
}


# Cat with no separator --------------------------------------------------------

cat0 <- function(...) {
  cat(..., sep = "")
}


# Capitalise the first letter of each string in a vector of strings ------------

firstCap <- function(str) {
  paste(toupper(substring(str, 1, 1)), substring(str, 2), sep = "")
}

# Return multiple objects from function as a list ------------------------------
## credit: G. Grothendieck
##
## Usage:
## fun <- function () {
##   ...
##   return(list(var1, var2, var3))
## }
## list[main1, main2, main3] <- fun()
## or
## list[main1, , main3] <- fun()

list <- structure(NA, class = "result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along = args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a = a,v = value[[i]])))
  }
  x
}

# Return only matching directories, not files ----------------------------------
## N.B. does not work for recursion!
## credit: Joshua Ulrich
my.list.dirs <- function(
  path = ".",         #
  pattern = NULL,     #
  all.dirs = FALSE,   #
  full.names = FALSE, #
  ignore.case = FALSE #
) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path = path, pattern = pattern, all.files = all.dirs,
                    full.names = TRUE, recursive = FALSE,
                    ignore.case = ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}


# Faster setdiff ---------------------------------------------------------------

my.setdiff <- function(x,y) {
  x[is.na(match(x,y))]
}


# Count number of lines of file ------------------------------------------------

my.count_lines <- function(path) {
  read.table(pipe(paste0("wc -l ", path))) [[1]]
}


# Mode -------------------------------------------------------------------------
## Returns first-appearing (not necessarily unique!) value of the set of modes

my.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# print large numbers as human-readable ----------------------------------------

my.human_readable <- function (x, digits = 1) {

  exponent <- seq(0, 9, by = 3)
  interval <- as.integer(10^exponent)
  suffix <- c("", "thousand", "million", "billion")

  i <- findInterval(x, interval)

  string <- paste(round(x/interval[i], digits), suffix[i])

  return(string)

}


# Grep multiple patterns at the same time --------------------------------------

## Verify if string contains all patterns, in any order and/or overlapping.
multi.grep <- function (
  pattern,        # string or vector of strings to be found in x
  x,              # string or vector of strings
  value = FALSE,  # return matching (default) or non-matching (invert == TRUE) elements?
  invert = FALSE, # return elements that do (FALSE) or do not (TRUE) match patterns
  ...             # arguments passed to grepl()
) {

  # compute TRUE or FALSE for each element in x
  res <- Reduce(`&`, lapply(pattern, function(pat) grepl(pattern = pat, x = x, ...)))

  # return matching (default) or non-matching (invert == TRUE) elements?
  if (invert)
    res <- !res

  # return match indices (default) or values (value == TRUE)
  res <- if (value) x[res] else which(res)

  return(res)

}
