#' Fit isotopic distribution-based segment-level model
#' 
#' @param observed_spectra xx
#' @param peptides_cluster xx
#' @param time_0_data xx
#' @param undeuterated_dists xx
#' @param starting_point xx
#' @param method xx
#' @param use_analytical_gradient xx
#' 
#' @export
#' 
fitIsoSegmentModel = function(observed_spectra,
                              peptides_cluster,
                              time_0_data,
                              undeuterated_dists,
                              starting_point = NULL,
                              method = "OLS",
                              use_analytical_gradient = TRUE,
                              max_iter = 100,
                              tolerance = 1e-2) {
  ols_optim_problem = getOptimizationProblem(observed_spectra, peptides_cluster, 
                                             time_0_data, undeuterated_dists,
                                             weights = NULL, theta = 1)
  num_parameters = getParametersCounts(peptides_cluster)
  starting_point = getStartingPoint(starting_point, num_parameters)
  peptide_segment_structure = getPeptideSegmentStructure(peptides_cluster)
  analytical_gradient = getAnalyticalGradient(use_analytical_gradient, observed_spectra,
                                       peptides_cluster, peptide_segment_structure,
                                       num_parameters, time_0_data, undeuterated_dists,
                                       weights = NULL, theta = 1)
  ols_solution = getOptimProblemSolution(starting_point, ols_optim_problem, analytical_gradient) # todo: handle non-convergence
  plgls_solution = getFinalSolution(method, use_analytical_gradient,
                                    starting_point, ols_solution,
                                    observed_spectra, 
                                    peptides_cluster,
                                    peptide_segment_structure, num_parameters,
                                    time_0_data,
                                    undeuterated_dists,
                                    max_iter,
                                    tolerance)
  plgls_solution
}

getParametersCounts = function(peptides_cluster) {
  segment_structure = unique(peptides_cluster[, list(Segment, MaxUptake)]) # change column name in the future
  segment_structure[["MaxUptake"]]
}

getStartingPoint = function(starting_point, num_parameters) {
  if (is.null(starting_point)) {
    starting_point = runif(num_parameters, -1e-3, -1e-5)
  }
  starting_point
}

getPeptideSegmentStructure = function(peptides_cluster) {
  segments = unique(peptides_cluster[, .(Segment, MaxUptake)])
  peptides_cluster[["Present"]] = 1
  ps_m = data.table::dcast(peptides_cluster, Peptide ~ Segment, value.var = "Present", fill = 0)
  ps_m = ps_m[, c("Peptide", as.character(segments$Segment)), with = FALSE]
  ps_m
}

getAnalyticalGradient = function(use_analytical_gradient, observed_spectra,
                                 peptides_cluster, peptide_segment_structure,
                                 num_parameters, time_0_data, undeuterated_dists,
                                 weights, theta) {
  if (use_analytical_gradient) {
    get_analytical_gradient_loss(observed_spectra, 
                                 peptides_cluster,
                                 peptide_segment_structure, num_parameters,
                                 time_0_data,
                                 undeuterated_dists,
                                 weights, theta)     
  } else {
    NULL
  } 
}

getOptimProblemSolution = function(starting_point, optim_problem, analytical_gradient) {
  solution = optim(starting_point, optim_problem, gr = analytical_gradient,
                   method = "BFGS")
  solution # todo: better processing
}

getFinalSolution = function(method, 
                            use_analytical_gradient,
                            starting_point, ols_solution,
                            observed_spectra, 
                            peptides_cluster,
                            peptide_segment_structure, num_parameters,
                            time_0_data,
                            undeuterated_dists,
                            max_iter,
                            tolerance) {
  if (method != "OLS") {
    initial_solution = starting_point
    current_solution = ols_solution$par
    current_parameters = current_solution
    
    iter = 1
    while(iter <= max_iter & max(abs(initial_solution - current_solution) / abs(initial_solution)) >= tolerance) {
      print(paste("Iteration: ", iter, ",", "relative difference:", max(abs(initial_solution - current_solution) / abs(initial_solution))))
      current_weights = getCurrentWeights(current_solution, observed_spectra, peptide_segment_structure, num_parameters, undeuterated_dists)
      current_theta = getCurrentTheta(current_solution, observed_spectra, current_weights, peptide_segment_structure, num_parameters, undeuterated_dists)
      
      current_optim_problem = getOptimizationProblem(observed_spectra, peptides_cluster, 
                                                     time_0_data, undeuterated_dists,
                                                     weights = current_weights, theta = current_theta)
      l2l = getAnalyticalGradient(use_analytical_gradient, observed_spectra,
                                  peptides_cluster, peptide_segment_structure,
                                  num_parameters, time_0_data, undeuterated_dists,
                                  current_weights, current_theta)
      
      current_optimized = optim(current_solution, current_optim_problem, 
                                gr = l2l,
                                method = "BFGS")
      initial_solution = current_solution
      current_solution = current_optimized$par
      iter = iter + 1
    }
    
    if (iter == 1) {
      current_optimized = ols_solution
      current_weights = NULL
      current_theta = 1
    }
    optimized = c(current_optimized, converged = !(max(abs(initial_solution - current_solution) / abs(initial_solution)) >= 1e-2), num_iter = iter - 1,
                  weights = current_weights, theta = current_theta)
    result = makeModelFittingOutput(observed_spectra, time_0_data, undeuterated_dists,
                                    peptide_segment_structure, num_parameters,
                                    optimized[["par"]], optimized[["theta"]], 
                                    optimized)
    result$OLS = ols_solution
    result
  } else {
    optimized = c(ols_solution, converged = ols_solution[["convergence"]] == 0, num_iter = 1,
                  weights = NULL, theta = NULL)
    list(FinalComparison = NULL,
         FittedProbabilities = NULL,
         OptimizationResult = optimized,
         OLS = ols_solution)
  }
}

makeModelFittingOutput = function(observed_spectra, time_0_data, undeuterated_dists,
                                  peptide_segment_structure,
                                  num_parameters, parameters, theta, optim_output) {
  peptides = peptide_segment_structure[["Peptide"]]
  time_points = unique(observed_spectra[["Time"]])
  
  seg_pars_final = getSegmentParametersFromBetas(parameters, num_parameters + 1)
  seg_probs_final = getSegmentProbabilitiesFromParams(seg_pars_final, time_points)
  pept_probs_final = getPeptideProbabilities(peptide_segment_structure, seg_probs_final)
  
  final_peak_ders = getFinalPeakDerivatives(observed_spectra,
                                            peptides_cluster,
                                            peptide_segment_structure,
                                            num_parameters,
                                            time_0_data,
                                            undeuterated_dists,
                                            parameters)
  
  final_theoretical_spectra = getExpectedSpectra(parameters, peptide_segment_structure, num_parameters + 1, observed_spectra, undeuterated_dists)
  final_comparison = merge(observed_spectra, final_theoretical_spectra,
                        by = c("Peptide", "Time", "IntDiff"),
                        sort = FALSE)
  
  probs_with_conf_ints = getProbabilitiesConfidenceIntervals(final_comparison, 
                                                             final_peak_ders, 
                                                             peptide_segment_structure, 
                                                             parameters, 
                                                             num_parameters, theta)
  list(FinalComparison = final_comparison,
       FittedProbabilities = probs_with_conf_ints,
       OptimizationResult = optim_output) # Todo: optimization history
}

getFinalPeakDerivatives = function(observed_spectra,
                                   peptides_cluster,
                                   peptide_segment_structure,
                                   num_parameters,
                                   time_0_data,
                                   undeuterated_dists,
                                   parameters) {
  peak_derivatives = getExpectedPeaksDerivatives(observed_spectra,
                                                 peptides_cluster,
                                                 peptide_segment_structure,
                                                 num_parameters,
                                                 time_0_data,
                                                 undeuterated_dists)
  final_peak_ders = peak_derivatives(parameters)
  final_peak_ders = final_peak_ders[!is.na(ExpectedPeak) & !is.na(Intensity) & ExpectedPeak > 0]
  final_peak_ders
}

#' @importFrom data.table rbindlist
getProbabilitiesConfidenceIntervals = function(final_comparison, 
                                               final_peak_ders, 
                                               peptide_segment_structure, 
                                               parameters, 
                                               num_parameters, theta) {
  time_points = unique(final_comparison[["Time"]])
  segments = colnames(peptide_segment_structure)[-1]
  parameters_covariance = getParametersCovarianceMatrix(final_comparison, final_peak_ders, parameters, theta) 
  probs_with_confidence = data.table::rbindlist(lapply(seq_along(time_points),
                                    function(ith_time) {
                                      getConfidenceIntervalSingleTimePoint(parameters_covariance, time_points, num_parameters,
                                                                           segments, parameters, nrow(final_comparison), ith_time)
                                    }))
  probs_with_confidence
}

getParametersCovarianceMatrix = function(final_comparison, final_peak_ders, parameters, theta) {
  p = length(parameters)
  final_comparison = final_comparison[!is.na(ExpectedPeak) & !is.na(Intensity) & ExpectedPeak > 0]
  var = final_comparison[, .(var = sum(((ExpectedPeak - Intensity) / (ExpectedPeak ^ theta))^2) / (.N - p))]$var
  
  par_ids = unique(final_peak_ders[, .(Segment, ParameterID)])
  par_ids[, ParamID := 1:nrow(par_ids)]
  
  ders = split(final_peak_ders, final_peak_ders[, .(Time, Peptide, IntDiff)])
  ders = ders[sapply(ders, function(x) !is.null(x) & nrow(x) > 0)]
  ders = lapply(ders, function(x) merge(x, par_ids, by = c("Segment", "ParameterID"),
                                        all.x = TRUE, all.y = TRUE))
  ders = lapply(ders, function(x) {
    x[, DerivativeEt := ifelse(is.na(DerivativeEt), 0, DerivativeEt)]
  })
  
  ders_matrices = lapply(ders, function(x) {
    x = x[, .(ParamID, DerivativeEt, ExpectedPeak)][order(ParamID)]
    (1 / (na.omit(unique(x$ExpectedPeak)))^(2*theta)) * (x$DerivativeEt %*% t(x$DerivativeEt))
  })
  
  total_matrix = matrix(0, nrow = p, ncol = p)
  for (m in ders_matrices) {
    total_matrix = total_matrix + m
  }
  cov_m = var * solve(total_matrix)
  cov_m
}


#' @importFrom numDeriv grad jacobian
getConfidenceIntervalSingleTimePoint = function(parameters_covariance, time_points, num_parameters,
                                                segments, parameters, n, time_id) {
  der_ith_time = getIthTimeDerivative(num_parameters, time_points, time_id)
  jac_time_ith = numDeriv::jacobian(der_ith_time, parameters)
  
  vars = diag(t(jac_time_ith) %*% parameters_covariance %*% jac_time_ith) / sqrt(n)
  
  seg_probs_ith = data.table(Segment = rep(segments, times = num_parameters + 1),
                            NumExchanged = unlist(lapply(num_parameters, function(x) 0:x)),
                            Probability = der_ith_time(parameters))
  
  seg_probs_ith[, Lower := Probability - qnorm(1 - 0.05/2) * sqrt(vars)]
  seg_probs_ith[, Upper := Probability + qnorm(1 - 0.05/2) * sqrt(vars)]
  seg_probs_ith[, Var := vars]
  seg_probs_ith[, Lower := max(0, Lower), by = c("Segment", "NumExchanged")]
  seg_probs_ith[, Upper := min(1, Upper), by = c("Segment", "NumExchanged")]
  seg_probs_ith[, Time := time_points[time_id]]
  seg_probs_ith
}

getIthTimeDerivative = function(num_parameters, times, i) {
  function(parameters) {
    segs = getSegmentParametersFromBetas(parameters, num_parameters + 1)
    probs = getSegmentProbabilitiesFromParams(segs, times)
    probs = probs[[i]]
    probs
    unlist(probs, F, F)
  }
}



