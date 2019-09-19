#make it work function for repeated try loop.----
makeitwork <- function(x, n = 4){
  for(k in 1:n){
    z <- try(x, silent = T)
    if(class(z) != 'try-error'){
      return(z)
      break
    }
    if(k == n){cat('All tries failed. bummer.\n')}
  }
}