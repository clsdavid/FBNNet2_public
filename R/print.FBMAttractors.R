#' A S3 method to print the FBM attractors
#' 
#' @param x The attactor details
#' @param ... other arguments.
#'@export
print.FBMAttractors <- function(x, ...) {
    print.FBMAttractors(x, ...)
}

# Custom print function for class FBMAttractor
print.FBMAttractors <- function(x, ...) {
    attractors <- x$Attractors
    genes <- x$Genes
    cat("Discovered Attractors via Fundamental Boolean Model ")
    cat(":\n")
    cat("Genes are encoded in the following order:")
    cat(":\n")
    cat(genes, sep = " ")
    cat(":\n\n")
    for (i in seq_along(attractors)) {
        if (length(attractors[[i]]) == 2) {
            statelen <- length(attractors[[i]][[1]])
            cat("Attractor ", i, " is a simple attractor consisting of ", length(attractors[[i]]) - 1, " state(s)", sep = "")
            cat(":\n\n")
            cat("|", "--<", rep("-", statelen), "|")
            cat("\n")
            cat("v", rep(" ", 2 + statelen), "^")
            cat("\n")
            cat(attractors[[i]][[1]], rep(" ", 3), "|")
            cat("\n")
            cat("|", rep(" ", 2 + statelen), "|")
            cat("\n")
            cat("v", rep(" ", 2 + statelen), "^")
            cat("\n")
            cat("|", rep("-", statelen), ">--", "|")
            cat("\n")
            cat("\n\n")
        } else {
            statelen <- length(attractors[[i]][[1]])
            cat("Attractor ", i, " is a complex attractor consisting of ", length(attractors[[i]]) - 1, " state(s)", sep = "")
            cat(":\n\n")
            cat("|", "--<", rep("-", statelen), "|")
            cat("\n")
            cat("v", rep(" ", 2 + statelen), "^")
            cat("\n")
            states <- attractors[[i]]
            j <- 1
            while (j < length(states)) {
                cat(states[[j]], rep(" ", 3), "|")
                cat("\n")
                cat("|", rep(" ", 2 + statelen), "|")
                cat("\n")
                j <- j + 1
            }
            cat("v", rep(" ", 2 + statelen), "^")
            cat("\n")
            cat("|", rep("-", statelen), ">--", "|")
            cat("\n")
            cat("\n\n")
        }
    }
    return(invisible(x))
}
