# tools.R - various utility functions

# function to add transparency to plot colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# transformations for logit-linear models 
logit <- function(x){
  y <- log(x/(1-x))
  return(y)
}
expit <- function(x){
  y <- exp(x)/(1+exp(x))
  return(y)
}
