
# Function created by T. Luke Smallman (t.l.smallman@ed.ac.uk)

ensemble_within_range<-function(target,proposal) {

   # Determine what proportion of a proposed PDF is within a target range
   # Returned value 0-1

   t_range = range(target, na.rm = TRUE)
   in_range = length(which(proposal >= t_range[1] & proposal <= t_range[2]))
   return(in_range / length(proposal))

} # ensemble_within_range
## Use byte compile
ensemble_within_range<-cmpfun(ensemble_within_range)
