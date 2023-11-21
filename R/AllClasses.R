# define S4 object class for scaDA

#' S4 class for storing scaDA dataset
#'
#' @description
#' This S4 class is designed to store the count matrix and associated column data information
#' for the scaDA package. currently the package can only applied to paired experiment design, single factor
#'
#' @slot data A matrix containing the read counts from experiments, must be positive integers
#' @slot colData A data.frame containing column metadata associated with experiment design.
setClass("scaDAdataset",
         slots = list(
           count = "matrix",  # Changed from 'data' to 'count'
           colData = "data.frame",
           params = "list"
         ),
         prototype = list(
           params = list()
         )
)


validity_scaDAdataset <- function(object) {
  if (!all(object@count >= 0) || !is.matrix(object@count) || !all(object@count == as.integer(object@count))) {
    return("Count slot must be a matrix of positive integers.")
  }

  if (ncol(object@count) != nrow(object@colData)) {
    return("Number of columns in count must match the number of rows in colData.")
  }

  TRUE
}

setValidity("scaDAdataset", validity_scaDAdataset)




#' scaDAdataset object and constructors; currently it need to be created from matrix
#' later we will make it seamlessly integerated with seurat package.


#' Constructor for scaDAdataset
#'
#' @param count A matrix containing the dataset.
#' @param colData A data.frame containing column metadata.
#' @return An object of class scaDAdataset.
#' @export
#' @import methods
scaDAdatasetFromMatrix <- function(count, colData) {
  new("scaDAdataset", count = count, colData = colData)
}

