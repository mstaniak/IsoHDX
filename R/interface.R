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
                              max_iter = 100) {
  ols_optim_problem = getOptimizationProblem(observed_spectra, peptides_cluster, 
                                             time_0_data, undeuterated_dists,
                                             weights = NULL, theta = 1)
  num_parameters = getParametersCounts(peptides_cluster)
  starting_point = getStartingPoint(starting_point, num_parameters)
  peptide_segment_structure = getPeptideSegmentStructure(peptides_cluster)
  analytical_gradient = getAnalyticalGradient(use_analytical_gradient, observed_spectra,
                                       peptides_cluster, peptide_segment_structure,
                                       num_parameters, time_0_data, undeuterated_dists,
                                       weights, theta)
  ols_solution = getOptimProblemSolution(starting_point, ols_optim_problem, analytical_gradient) # todo: handle non-convergence
  plgls_solution = getFinalSolution(method, use_analytical_gradient,
                                    starting_point, ols_solution,
                                    observed_spectra, 
                                    peptides_cluster,
                                    peptide_segment_structure, num_parameters,
                                    time_0_data,
                                    undeuterated_dists,
                                    max_iter)
  plgls_solution
}

getParametersCounts = function(peptides_cluster) {
  segment_structure = unique(peptides_cluster[, list(Segment, MaxUptake)]) # change column name in the future
  segment_structure[["MaxUptake"]]
}

getStartingPoint = function(starting_point, num_parameters) {
  if (is.null(starting_point)) {
    starting_point = runif(num_parameters, 1e-5, 1e-3)
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
                            max_iter) {
  if (method != "OLS") {
    initial_solution = starting_point
    current_solution = ols_solution$par
    current_parameters = current_solution
    
    iter = 1
    while(iter <= max_iter & sum((initial_solution - current_solution)^2) >= 1e-6) {
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
    c(current_optimized, converged = sum((initial_solution - current_solution)^2) >= 1e-6, num_iter = iter - 1)
  } else {
    c(ols_solution, converged = ols_solution[["convergence"]] == 0, num_iter = 1)
  }
}

# seg_pars_final = getSegmentParametersFromBetas(current_optimized$par, num_parameters + 1)
# seg_probs_final = getSegmentProbabilitiesFromParams(seg_pars_final, unique(observed_spectra$Time))
# pept_probs_final = getPeptideProbabilities(peptide_segment_structure, seg_probs_final)
# 
# # Gradient
# 
# peak_derivatives = get_expected_peaks_derivatives(observed_spectra, 
#                                                   peptides_cluster,
#                                                   peptide_segment_structure, num_parameters,
#                                                   peaks_0_7_info[Charge == 2],
#                                                   und_d)
# final_peak_ders = peak_derivatives(current_solution)
# final_peak_ders = final_peak_ders[!is.na(ExpectedPeak) & !is.na(Intensity) & ExpectedPeak > 0]
# 
# final_theoretical_spectra = getExpectedSpectra2(current_solution, peptide_segment_structure, num_parameters + 1, observed_spectra, undeuterated_dists)
# final_compare = merge(observed_spectra, final_theoretical_spectra,
#                       by = c("Peptide", "Time", "IntDiff"),
#                       sort = FALSE)
# theta = current_theta
# p = length(current_solution)
# final_compare = final_compare[!is.na(ExpectedPeak) & !is.na(Intensity) & ExpectedPeak > 0]
# var2 = final_compare[, .(var = sum(((ExpectedPeak - Intensity) / (ExpectedPeak ^ theta))^2) / (.N - p))]$var
# 
# der_check = final_peak_ders[, .(Time, Peptide, IntDiff, Segment, ParameterID, DerivativeEt, 
#                                 ExpectedPeak, Intensity)]
# par_ids = unique(der_check[, .(Segment, ParameterID)])
# par_ids[, ParamID := 1:nrow(par_ids)]
# 
# 
# ders = split(final_peak_ders, final_peak_ders[, .(Time, Peptide, IntDiff)])
# ders = ders[sapply(ders, function(x) !is.null(x) & nrow(x) > 0)]
# ders = lapply(ders, function(x) merge(x, par_ids, by = c("Segment", "ParameterID"),
#                                       all.x = T, all.y = T))
# ders = lapply(ders, function(x) {
#   x[, DerivativeEt := ifelse(is.na(DerivativeEt), 0, DerivativeEt)]
# })
# ders_matrices = lapply(ders, function(x) {
#   x = x[, .(ParamID, DerivativeEt, ExpectedPeak)][order(ParamID)]
#   (1 / (na.omit(unique(x$ExpectedPeak)))^(2*theta)) * (x$DerivativeEt %*% t(x$DerivativeEt))
# })
# 
# 
# get_ith_time_der = function(num_parameters, times, i) {
#   function(current_solution) {
#     segs = getSegmentParametersFromBetas(current_solution, num_parameters + 1)  
#     probs = getSegmentProbabilitiesFromParams(segs, times)
#     probs = probs[[i]]
#     probs
#     unlist(probs, F, F)
#   }
# }
# 
# 
# total_matrix = matrix(0, nrow = 18, ncol = 18)
# for (m in ders_matrices) {
#   total_matrix = total_matrix + m
# }
# 
# round(var2 * diag(solve(total_matrix)), 2)
# 
# cov_m = var2 * solve(total_matrix)
# det(cov_m)
# pracma::isposdef(cov_m)
# 
# der_2_time = get_ith_time_der(num_parameters, times, 2)
# der_2_time(current_solution)
# length(der_2_time(current_solution))
# 
# grad_time_2 = numDeriv::grad(der_2_time, current_solution)
# grad_time_2
# 
# jac_time_2 = numDeriv::jacobian(der_2_time, current_solution)
# 
# t(grad_time_2) %*% cov_m %*% t(t(grad_time_2))
# t(grad_time_2) %*% cov_m %*% grad_time_2
# 
# round(diag(t(jac_time_2) %*% cov_m %*% jac_time_2) / sqrt(nrow(final_compare)), 6)
# 
# num_parameters
# der_2_time(current_solution)
# vars = diag(t(jac_time_2) %*% cov_m %*% jac_time_2) / sqrt(nrow(final_compare))
# seg_probs_t2 = data.table(Segment = rep(as.character(segments$Segment), times = num_parameters + 1),
#                           NumExchanged = unlist(sapply(num_parameters, function(x) 0:x)),
#                           Probability = der_2_time(current_solution))
# seg_probs_t2[, Lower := Probability - qnorm(1 - 0.05/2) * sqrt(vars)]
# seg_probs_t2[, Upper := Probability + qnorm(1 - 0.05/2) * sqrt(vars)]
# seg_probs_t2[, Var := vars]
# seg_probs_t2[, Lower := max(0, Lower), by = c("Segment", "NumExchanged")]
# seg_probs_t2[, Upper := min(1, Upper), by = c("Segment", "NumExchanged")]
#
