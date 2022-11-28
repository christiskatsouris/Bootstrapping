#################################################################################################
# permutation_test.R
#
# A function for performing permutation testing that ensures that we sample
# from the space of all permutations without replacement.
#
# Reference: https://gist.github.com/Sarkom 
#################################################################################################

library(hash)
library(gmp)
library(digest)
library(sn)  # For producing skewed distributions


##########################

permutation.test <- function(sample1, sample2, metric, I = 1000)
{# begin-of-function
 
  # Perform a permutation test by sampling without replacement from the space of permutations
  start = Sys.time()
  
  # Get test statistic
  Tobs = metric(sample1) - metric(sample2)
  
  # Pool samples
  n1 = length(sample1)
  n2 = length(sample2)
  n = n1 + n2
  allsamples = c(sample1, sample2)
  
  # Do a bunch of iterations
  Tcalc = rep(0,I)
  combins = hash()
  
  for( i in 1:I )
  {
    # Generate random combinations
    choicevector = c(rep(0,n1),rep(1,n2))
    samplechoice = sample(choicevector,n)
    # Convert the samplechoice vector of 1s and 0s to a large binary number
    key = as.bigz(paste("0b",paste(samplechoice,collapse=""),sep=""))
    # Check if it is in the hash table - key length is max of 256 for the hash table
   
   if( n < 1200 )  # 256 = log_32(2^1280)
    {
      key = as.character(key, b = 32)
      while( !is.null( combins[[ key ]] ) ) 
      {
        cat("Collision!  Resampling...\n")
        samplechoice = sample(choicevector,n)
        key = paste(samplechoice,collapse="")    
        key = as.character(key, b = 32)
      }
      combins[[ key ]] = 1
    }
    else
    {
      key = digest(as.character(key,b=32), algo = "sha1")
      while( !is.null( combins[[ key ]] ) ) {
        cat("Collision!  Resampling...\n")        
        samplechoice = sample(choicevector,n)
        key = paste(samplechoice,collapse="")		
        key = digest(as.character(key,b=32), algo = "sha1")
      }
      combins[[ key ]] = 1
    }
    
    newsample1 = allsamples[samplechoice == 0]
    newsample2 = allsamples[samplechoice == 1]
    if  ((length(newsample1) != length(sample1)) |(length(newsample2) != length(sample2)) )
   {
      stop("Error in resampling!")
    }
    Tcalc[i] = metric(newsample1) - metric(newsample2)
  }
  
  # Determine p-value
  p = sum(abs(Tcalc) >= abs(Tobs)) / I
  
  end = Sys.time()
  # print(difftime(end, start))
  
  return(list(pvalue = p,
              p.value = p,
              Tobs = Tobs,
              numiter = I,
              time = difftime(end, start),
              Tcalc = Tcalc))
}# end-of-function

#################################################################################################
#################################################################################################

# Test 1: N(0,1) vs N(100,1)

x = rnorm(1000, 0, 1)
y = rnorm(1000, 100, 1)

cat("N(0,1) vs N(100,1) -- compare normal distributions with same variance and different means:\n")
cat("With var:")
replicate(5, permutation.test(x, y, var)$pvalue)
cat("With sd:")

replicate(5, permutation.test(x, y, sd)$pvalue)
