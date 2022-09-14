# Tools for simulating isotopic distributions of peptides from time-varying HDX-MS experiments

#' Simulate spectra for given sequences and time points
#' 
#' @param peptide_sequences character vector of peptide sequences
#' @param nums_exchangeable numeric vector of exchangeable hydrogens counts for each peptide
#' @param time_points numeric vector of time points
#' @param lambdas numeric vector of parameters for atom-specific exchange rates
#' @param total_intensity total intensities in a spectrum
#' @param error_sd error for intensities (standard deviation of N(0, error_sd))
#' 
#' @import data.table
#' 
#' @return data.table
#' 
getSpectra = function(peptide_sequences, nums_exchangeable,
                      time_points, lambdas,
                      total_intensity = 1e6, error_sd = 0.1) {
  undeuterated_distributions_raw = lapply(
    peptide_sequences, function() BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(x))
  )
  undeuterated_dists = lapply(undeuterated_dists_raw,
                              function(x) x$isoDistr)
  
  data.table::rbindlist(
    lapply(seq_along(undeuterated_distributions),
           function(i) {
             num_exchangeable = nums_exchangeable[i]
             monoisotopic_mass = undeuterated_distributions_raw[[i]]$monoisotopicMass
             undeuterated_iso = undeuterated_distributions[[i]]
             lambda = lambdas[[i]]
             
             result_matrix = matrix(c(undeuterated_iso,
                                      rep(0, num_exchangeable + 1)),
                                    nrow = num_exchangeable + 1,
                                    ncol = length(undeuterated_iso) + num_exchangeable,
                                    byrow=TRUE)
             
             simulated_intensities = genData_hetero_new(refM = result_matrix, 
                                                        totIntens = total_intensity, 
                                                        lambda = lambda,
                                                        timepoints= time_points, 
                                                        error = error_sd) # generates peaks
             
             times = rep(time_points, 
                         each = length(undeuterated_iso) + num_exchangeable)
             num_time_points = length(time_points)
             num_probabilities = dim(result_matrix)[2]
             
             mass = rep(monoisotopic_mass:(monoisotopic_mass + dim(result_matrix)[2] - 1),
                        times = num_time_points)
             
             data.table::data.table(Peptide = peptide_sequences[i],
                                    Time = times,
                                    Mass = mass,
                                    Intensity = simulated_intensities)
           }))
}


#' @keywords internal
#' @author Jurgen Claesen
genData_hetero_new <- function(refM, totIntens, lambda, timepoints, error){
  
  intens <- c()
  E <- rep(0, dim(refM)[2])
  for(i in 1:length(timepoints)){
    tau <- timepoints[i]
    intm <- calcProbabilities(lambda, tau)
    intm2 <- totIntens*(intm%*%refM)
    E <- error*intm2*rnorm(n=dim(refM)[2],mean=0,sd=1)
    intens <- c(intens,intm2+E)
  }
  intens
}

#' @keywords internal
#' @author Jurgen Claesen
calcProbabilities = function(lambda, tau){
  tmp <- exp(-lambda*tau)
  exp_theta <- tmp/(1-tmp)
  
  n <- length(lambda)
  delta <- rep(0,n+1)
  delta[n+1] <- prod(1-tmp)
  coeff <- getElementPolynom(exp_theta,n)
  delta[1:n] <- coeff[n:1]*delta[n+1]
  return(delta)
}

#' @keywords internal
getElementPolynom <- function(exp_theta,n){
  elem <- rep(0,n)
  elem[1] <- sum(exp_theta)
  p <- rep(0,n)
  
  for(j in 1:n){
    p[j] <- sum(exp_theta^j)
  }
  
  for(k in 2:n){
    tmp <- rep(0,k-1)
    for(i in 1:(k-1)){
      tmp[i] <- (-1)^(i+1)*elem[k-i]*p[i]
    }
    elem[k] <- 1/k*(sum(tmp)+(-1)^(k+1)*p[k])
  }
  return(elem)
}
