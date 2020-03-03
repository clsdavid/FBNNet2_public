#'@export
print.FundamentalBooleanNetwork <- function(x, ...) {
    print.FundamentalBooleanNetwork(x, ...)
}

# Custom print function for class FundamentalBooleanNetwork
#'@export
print.FundamentalBooleanNetwork <- function(x, ...) {
    cat("Fundamental Boolean Network with ", length(x$genes), "genes\n")
    cat("Genes involved:\n", paste(x$genes, collapse = ", ", sep = ", "), "\n", sep = "")
    cat("\nNetworks:")
    mapply(function(gene, interaction) {
        cat("\nMultiple Transition Functions for ", gene, " with decay value = ", x$timedecay[[gene]], ":\n", sep = "")
        # print original expressions read from the files (if available)
        lapply(names(interaction), function(funcName) {
            cat(funcName, ": ", gene, " = ", interaction[[funcName]]$expression, sep = "")
            
            if (!is.null(interaction[[funcName]]$error)) {
                cat(" (")
                
                # if (!is.null(interaction[[funcName]]$error)) cat('Error: ',as.character(as.numeric(interaction[[funcName]]$error) * 100),
                # '%',sep='')
                
                # cat(', ')
                
                if (!length(interaction[[funcName]]$probability) == 0 | !is.null(interaction[[funcName]]$probability)) 
                  cat("Confidence: ", as.character(interaction[[funcName]]$probability), sep = "")
                
                cat(", ")
                
                if (!length(interaction[[funcName]]$timestep) == 0 | !is.null(interaction[[funcName]]$timestep)) 
                  cat("TimeStep: ", interaction[[funcName]]$timestep, sep = "") else cat("TimeStep: ", 1, sep = "")
                
                cat(")")
            }
            cat("\n")
        })
    }, names(x$interactions), x$interactions)
    
    if (sum(x$fixed != -1) > 0) {
        cat("\nKnocked-out and over-expressed genes:\n")
        mapply(function(gene, fixed) {
            if (fixed != -1) 
                cat(gene, " = ", fixed, "\n", sep = "")
        }, x$genes, x$fixed)
    }
    
    return(invisible(x))
}
