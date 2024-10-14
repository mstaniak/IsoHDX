#' Create optimization problem to estimate segment-specific exchange rates via OLS
#' 
#' @param observed_spectra data.table with observed isotopic distributions
#' @param peptides_cluster data.table with information about peptides cluster
#' @param time_0_data data.table with information about isotopic distributions of 
#' undeuterated peptides
#' @param weights optional data.table with weights for each time point, peptide and isotopic peak
#' @param theta optional scalar variance parameter (see more in the Details section)
#' 
#' @import data.table
#' @export
#' 
getOptimizationProblem = function(observed_spectra,
                                  peptides_cluster,
                                  time_0_data, undeuterated_dists,
                                  weights = NULL, theta = 1) {
  times = unique(observed_spectra$time)
  segments = unique(peptides_cluster[, .(Segment, MaxUptake)])
  peptides_cluster[, Present := 1]
  ps_m = data.table::dcast(peptides_cluster, Peptide ~ Segment, value.var = "Present", fill = 0)
  ps_m = ps_m[, c("Peptide", as.character(segments$Segment)), with = FALSE]
  num_parameters = segments[["MaxUptake"]]
  
  function(parameters) {
    theoretical_spectra = getExpectedSpectra(parameters,
                                             ps_m,
                                             num_parameters + 1,
                                             observed_spectra,
                                             undeuterated_dists)
    
    compare = merge(theoretical_spectra,
                    observed_spectra,
                    by = c("Peptide", "Time", "IntDiff"),
                    all.x = T, all.y = T)
    # compare[, ExpectedPeak := ifelse(is.na(ExpectedPeak), 0, ExpectedPeak)] # based on comments in the draft
    # compare[, Intensity := ifelse(is.na(Intensity), 0, Intensity)] # based on comments in the draft
    if (!is.null(weights)) {
      compare = merge(compare, weights,
                      by = c("Peptide", "Time", "IntDiff"),
                      # by = c("Peptide", "Charge", "Time", "IntDiff"),
                      all.x = T, all.y = T)
      sum( ((compare[["Weight"]] ^ theta) * (compare[["Intensity"]] - compare[["ExpectedPeak"]])) ^ 2,
           na.rm = TRUE)
    } else {
      sum( ((compare[["Intensity"]] - compare[["ExpectedPeak"]])) ^ 2,
           na.rm = TRUE)
    }
  }
}

#' @keywords internal
getSegmentParametersFromBetas = function(betas, num_parameters) {
  all_probs = lapply(seq_along(num_parameters), function(i) {
    if (i == 1) {
      subset_betas = betas[1:num_parameters[i]]
    } else {
      subset_betas = betas[(sum(num_parameters[1:(i - 1)]) + 1):(sum(num_parameters[1:i]))]
    }
    subset_betas
  })
  all_probs
}

#' @keywords internal
getSegmentProbabilitiesFromParams = function(segment_params, time_points) {
  lapply(time_points, function(time) {
    lapply(segment_params, getProbabilitiesFromBetas, time = time)
  })
}

#' @keywords internal
getExpectedPeakHeights = function(total,
                                  probabilities,
                                  undeuterated_dist,
                                  num_exchangeable) {
  
  sapply(1:(length(undeuterated_dist) + num_exchangeable), 
         function(k) {
           sum(sapply(max(1, k - num_exchangeable):min(length(undeuterated_dist), k),
                      function(d) {
                        probabilities[k - d + 1] * undeuterated_dist[d] * total
                      }), na.rm = TRUE)
         })
}

#' @keywords internal
getPeptideProbabilities = function(pept_seg_struct, segment_probs) {
  lapply(segment_probs, function(probs_in_time) {
    lapply(seq_len(nrow(pept_seg_struct)), function(ith_row) {
      segments = as.logical(unlist(pept_seg_struct[ith_row, colnames(pept_seg_struct) != "Peptide", with = FALSE]))
      pept_segments = probs_in_time[segments]
      num_exchangeable = sum(sapply(pept_segments, length)) - length(pept_segments)
      getExchangeProbabilities(pept_segments, num_exchangeable)
    })
  })
}

#' @keywords internal
getExpectedSpectra = function(parameters,
                              pept_seg_struct,
                              param_counts,
                              observed_spectra,
                              undeuterated_dists) {
  times = unique(observed_spectra$Time)
  params_by_seq = getSegmentParametersFromBetas(parameters, param_counts)
  probs_by_time_seg = getSegmentProbabilitiesFromParams(params_by_seq, times)
  probs_by_time_seg = lapply(probs_by_time_seg, function(x) lapply(x, function(y) {
    y = ifelse(is.nan(y), (1 - sum(y[!is.nan(y)])) / sum(is.nan(y)), y)
    y = y / sum(y)
    y
  }))
  pept_probs = getPeptideProbabilities(pept_seg_struct, probs_by_time_seg)
  
  rbindlist(lapply(seq_along(pept_probs), function(ith_time) {
    probs_in_time = pept_probs[[ith_time]]
    rbindlist(lapply(seq_along(probs_in_time), function(ith_peptide) {
      probs = probs_in_time[[ith_peptide]]
      undeuterated_probs = undeuterated_dists[[unique(pept_seg_struct$Peptide[ith_peptide])]]
      peaks_heights = getExpectedPeakHeights(1, probs$Probability, undeuterated_probs, max(probs$NumExchanged))
      list(Peptide = unique(pept_seg_struct$Peptide[ith_peptide]),
           Time = times[ith_time],
           IntDiff = 0:(length(peaks_heights) - 1),
           ExpectedPeak = peaks_heights)
    }))
  }))
}

#' @keywords internal
getProbabilitiesFromBetas = function(betas, time) {
  intercept = betas[1]
  betas = betas[-1]
  # betas = betas - max(betas) 
  total = sum(exp(betas * time + intercept)) + 1
  probs = exp(betas * time + intercept) / total
  probs = c(probs, 1 - sum(probs))
  if (any(is.nan(probs))) browser()
  probs
}

#' @keywords internal
getProtonMass = function() {
  1.0072826748
}

#' @keywords internal
get_analytical_derivatives = function(peptides_cluster, pept_seg_struct, undeuterated_dists, num_parameters, time_points, parameters) {
  total_num_params = sum(num_parameters)
  n_pars = num_parameters
  names(n_pars) = colnames(pept_seg_struct)[-1]
  all_segments_names = names(n_pars)
  
  result = lapply(seq_len(nrow(pept_seg_struct)), function(peptide_id) {
    peptide = pept_seg_struct[["Peptide"]][peptide_id]
    segments_in_peptide = as.logical(unlist(pept_seg_struct[peptide_id, -1], F, F))
    segment_params = getSegmentParametersFromBetas(parameters, num_parameters + 1)
    segment_probabilities = getSegmentProbabilitiesFromParams(segment_params, time_points)
    n = sum(n_pars[all_segments_names][segments_in_peptide])
    
    lapply(seq_along(time_points), function(time_id) {
      time = time_points[time_id]
      probs = segment_probabilities[[time_id]]
      peptide_probs = getExchangeProbabilities(probs[segments_in_peptide],
                                               sum(n_pars[all_segments_names][segments_in_peptide]),
                                               FALSE)
      
      lapply(peptide_probs$NumExchanged, function(j) {
        lapply(seq_along(segment_params[segments_in_peptide]), function(segment_id) {
          present_segments = lapply(seq_along(segments_in_peptide), function(segment_id) unique(segments_in_peptide[segment_id]))
          segment_to_remove = which(segments_in_peptide)[segment_id]
          present_segments[[segment_to_remove]] = FALSE
          present_segments = unlist(present_segments, F, F)
          
          if (any(present_segments)) {
            segment_probabilities_other = probs[present_segments]
            
            num_exchangeable = sum(num_parameters[present_segments])
            num_exch_segment_k = sum(num_parameters[segments_in_peptide]) - num_exchangeable
            other_segments_prob = getExchangeProbabilities(segment_probabilities_other,
                                                           num_exchangeable,
                                                           FALSE)
            num_parameters_segment = length(segment_params[segments_in_peptide][[segment_id]])
            
            min_id = max(c(0, j - (n - num_exch_segment_k)))
            max_id = min(c(j, num_exch_segment_k))
            
            rbindlist(lapply(seq_len(num_parameters_segment), function(l) {
              list(Peptide = peptide,
                   NumExchanged = j,
                   Time = time,
                   Segment = colnames(pept_seg_struct[, -1, with = FALSE])[segments_in_peptide][segment_id],
                   ParameterID = l,
                   ForDerivative = sum(sapply(min_id:max_id, function(i) {
                     get_derivative(segment_params[segments_in_peptide][[segment_id]], time, 
                                    l, i, num_parameters_segment - 1) * other_segments_prob$Probability[other_segments_prob$NumExchanged == j - i]
                   })))
            }))            
          } else {
            num_exchangeable = sum(num_parameters[present_segments])
            num_exch_segment_k = sum(num_parameters[segments_in_peptide]) - num_exchangeable
            num_parameters_segment = length(segment_params[segments_in_peptide][[segment_id]])
            
            min_id = max(c(0, j - (n - num_exch_segment_k)))
            max_id = min(c(j, num_exch_segment_k))
            
            rbindlist(lapply(seq_len(num_parameters_segment), function(l) {
              list(Peptide = peptide,
                   NumExchanged = j,
                   Time = time,
                   Segment = colnames(pept_seg_struct[, -1, with = FALSE])[segments_in_peptide][segment_id],
                   ParameterID = l,
                   ForDerivative = sum(sapply(min_id:max_id, function(i) {
                     get_derivative(segment_params[segments_in_peptide][[segment_id]], time, 
                                    l, i, num_parameters_segment - 1) * 1
                   })))
            }))   
          }
        }) 
      })
    })
  })
  rbindlist(unlist(unlist(unlist(result, F, F), F, F), F, F))
}

#' @keywords internal
getUndeuteratedDists = function(time_0_data) {
  undeuterated_dists = lapply(seq_len(nrow(time_0_data)), function(i) {
    x =  BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(time_0_data$Sequence[i]),
                         nrPeaks = time_0_data$NumPeaks[i])$isoDistr
    x
  })
  names(undeuterated_dists) = time_0_data$Sequence
  undeuterated_dists
}

#' Analytical gradient of PLS-GLS function
#' 
#' @param observed_spectra description
#' @param peptides_cluster description
#' @param pept_seg_struct description
#' @param num_parameters description
#' @param time_0_data description
#' @param undeuterated_dists description
#' @param weights description
#' @param theta description
#' 
#' @export
#' 
get_analytical_gradient_loss = function(observed_spectra, peptides_cluster,
                                        pept_seg_struct, num_parameters,
                                        time_0_data, undeuterated_dists,
                                        weights = NULL, theta = 1) {
  time_points = unique(observed_spectra$Time)
  observed_spectra = observed_spectra
  peptides_cluster = peptides_cluster
  pept_seg_struct = pept_seg_struct
  num_parameters = num_parameters
  time_0_data = time_0_data
  undeuterated_dists = undeuterated_dists
  num_ex_peptide = peptides_cluster[, .(Peptide, Segment, MaxUptake)]
  num_ex_peptide = num_ex_peptide[, .(NumEx = sum(MaxUptake)), by = "Peptide"]
  
  function(parameters) {
    analytical_gradient_probs = get_analytical_derivatives(peptides_cluster, pept_seg_struct, undeuterated_dists, num_parameters, time_points, parameters)
    theoretical_spectra = getExpectedSpectra(parameters, pept_seg_struct, num_parameters + 1, observed_spectra, undeuterated_dists)
    compare = merge(observed_spectra, theoretical_spectra,
                    by = c("Peptide", "Time", "IntDiff"),
                    sort = FALSE)
    
    derivatives = lapply(unique(analytical_gradient_probs$Time), function(time) {
      lapply(unique(analytical_gradient_probs$Segment[analytical_gradient_probs$Time == time]), function(segment) {
        lapply(unique(analytical_gradient_probs$ParameterID[analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment]), function(param_id) {
          lapply(unique(analytical_gradient_probs$Peptide[analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment]), function(peptide) {
            undeuterated_probabilities_current = undeuterated_dists[[peptide]]
            num_exchangeable = num_ex_peptide$NumEx[num_ex_peptide$Peptide == peptide]
            
            lapply(observed_spectra[observed_spectra$Peptide == peptide & observed_spectra$Time == time, unique(Charge)], function(charge) {
              total = sum(observed_spectra$Intensity[observed_spectra$Peptide == peptide & observed_spectra$Time == time & observed_spectra$Charge == charge], na.rm = T)
              lapply(1:(length(undeuterated_probabilities_current) + num_exchangeable), 
                     function(k) {
                       peak_range = max(1, k - num_exchangeable):min(length(undeuterated_probabilities_current), k)
                       range_len = length(peak_range)
                       list(Peptide = peptide,
                            Segment = segment,
                            Time = time,
                            Charge = charge,
                            ParameterID = param_id,
                            DerivativeEt = sum(vapply(peak_range,
                                                      function(d) {
                                                        filt = analytical_gradient_probs$ParameterID == param_id & analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment & analytical_gradient_probs$Peptide == peptide & analytical_gradient_probs$NumExchanged == k - d
                                                        x = undeuterated_probabilities_current[d] * total * analytical_gradient_probs$ForDerivative[filt]
                                                        x
                                                      }, numeric(1)), na.rm = TRUE),
                            IntDiff = k - 1)
                     })
            })
          })
        })
      })
    })
    derivatives_df = rbindlist(unlist(unlist(unlist(unlist(unlist(derivatives, F, F), F, F), F, F), F, F), F, F))
    
    res = merge(derivatives_df, compare, 
                by = c("Time", "Peptide", "Charge", "IntDiff"),
                allow.cartesian = T, sort = FALSE)
    
    if (!is.null(weights)) {
      res = merge(res, weights,
                  by = c("Peptide", "Time", "IntDiff"), sort = FALSE)
      res = res[, .(Derivative = 2 * sum((Weight^(2*theta)) * (ExpectedPeak - Intensity) * DerivativeEt, na.rm = T)), by = c("Segment", "ParameterID")]
    } else {
      res = res[, .(Derivative = 2 * sum((ExpectedPeak - Intensity) * DerivativeEt, na.rm = T)), by = c("Segment", "ParameterID")]
    }
    res[, Derivative]
  }
}

#' @keywords internal
getExchangeProbabilities = function(by_segment_probabilities,
                                    num_exchangeable,
                                    approximate_root = FALSE) {
  first_nonzero = vapply(by_segment_probabilities, function(x) min(which(x >= 1e-8)), numeric(1))
  last_nonzero = vapply(by_segment_probabilities, function(x) max(which(abs(x) >= 1e-8)), numeric(1))
  lengths = vapply(by_segment_probabilities, length, numeric(1))
  nonzero_start = sum(first_nonzero) - length(first_nonzero)
  nonzero_end = sum(last_nonzero != lengths)
  by_segment_probabilities = lapply(seq_along(by_segment_probabilities), function(i) by_segment_probabilities[[i]][first_nonzero[i]:last_nonzero[i]])
  no_exchange_probability = prod(sapply(by_segment_probabilities, function(x) x[1]))
  num_exchangeable_complete = num_exchangeable
  num_exchangeable = num_exchangeable_complete - nonzero_start - nonzero_end
  power_sums = lapply(by_segment_probabilities, function(x) {
    power_sums = vector("numeric", num_exchangeable + 1)
    power_sums[1] = no_exchange_probability
    coefs = c(x, rep(0, num_exchangeable + 1 - length(x)))
    
    for (i in 2:(num_exchangeable + 1)) {
      if (i <= length(x)) {
        if (i == 2) {
          power_sums[i] = -(i - 1) * coefs[i] / coefs[i - 1] 
        } else {
          power_sums[i] = ((-(i - 1) * coefs[i]) - sum(rev(coefs[2:(i - 1)]) * power_sums[2:(i - 1)])) / coefs[1]
        }
      } else {
        power_sums[i] = -sum(rev(coefs[2:(i - 1)]) * power_sums[2:(i - 1)]) / coefs[1]
      }
    }
    power_sums
  })
  roots_list = colSums(matrix(unlist(power_sums), nrow = length(power_sums), byrow = T))[-1]
  
  polynomials = vector("list", num_exchangeable)
  for (i in seq_len(num_exchangeable)) {
    if (i == 1) {
      polynomials[[i]] = - no_exchange_probability * roots_list[1]
    } else {
      polynomials[[i]] = (-1 / i) * sum(rev(c(no_exchange_probability, unlist(polynomials[1:(i - 1)], FALSE, FALSE))) * roots_list[1:i])
    }
  }
  probabilities = lapply(polynomials, Re)
  probabilities = c(no_exchange_probability, probabilities)
  probabilities = unlist(probabilities, FALSE, FALSE)
  
  list(NumExchanged = 0:num_exchangeable_complete,
       Probability = c(rep(0, nonzero_start), probabilities, rep(0, nonzero_end)))
}

#' @keywords internal
get_derivative = function(betas_segment, time, l, j, n_k) {
  index_r = l
  l = l - 2
  if (l == -1 & j < n_k) {
    intercept = betas_segment[1]
    betas = betas_segment[-1]
    total = 1 + sum(exp(betas * (time - intercept)))
    current_beta = betas_segment[index_r]
    this_exp = exp(current_beta * (time - intercept))
    jth_exp = exp(betas[j + 1] * (time - intercept))
    jth_exp * (-betas[j + 1] + sum((betas - betas[j + 1]) * exp(betas * (time - intercept)))) / total ^ 2
  } else if (l == -1 & j == n_k) {
    intercept = betas_segment[1]
    betas = betas_segment[-1]
    total = 1 + sum(exp(betas * (time - intercept)))
    numerator = sum(-betas * exp(betas * (time - intercept)))
    - numerator / total ^ 2
  } else if (l > -1 & l == j & j < n_k) {
    intercept = betas_segment[1]
    betas = betas_segment[-1]
    total = 1 + sum(exp(betas * (time - intercept)))
    current_beta = betas_segment[index_r]
    this_exp = exp(current_beta * (time - intercept))
    (time - intercept) * this_exp * (total - this_exp) / total ^ 2
  } else if (l > -1 & l != j & j < n_k) {
    intercept = betas_segment[1]
    betas = betas_segment[-1]
    total = 1 + sum(exp(betas * (time - intercept)))
    current_beta = betas_segment[index_r]
    this_exp = exp(current_beta * (time - intercept))
    jth_exp = exp(betas[j + 1] * (time - intercept))
    -(time - intercept) * jth_exp * this_exp / total ^ 2
  } else if (l > -1 & l < j & j == n_k) {
    intercept = betas_segment[1]
    betas = betas_segment[-1]
    total = 1 + sum(exp(betas * (time - intercept)))
    current_beta = betas_segment[index_r]
    this_exp = exp(current_beta * (time - intercept))
    -(time - intercept) * this_exp / total ^ 2
  }
}

#' @keywords internal
getCurrentWeights = function(current_parameters, observed_spectra, pept_seg_struct, num_parameters, undeuterated_dists) {
  expected_spectra = getExpectedSpectra(current_parameters, pept_seg_struct,
                                        num_parameters + 1, observed_spectra, undeuterated_dists)
  expected_spectra = expected_spectra[ExpectedPeak > 0 & !is.na(ExpectedPeak)]
  prod = (1 / nrow(expected_spectra)) * sum(log(expected_spectra$ExpectedPeak))
  grand_mean = exp(prod)
  expected_spectra[, Weight := grand_mean / (ExpectedPeak)]
  expected_spectra[, ExpectedPeak := NULL]
  expected_spectra
}

#' @keywords internal
getCurrentTheta = function(current_parameters, observed_spectra, weights, pept_seg_struct, num_parameters, undeuterated_dists) {
  theoretical_spectra = getExpectedSpectra(current_parameters, pept_seg_struct,
                                           num_parameters + 1, observed_spectra, undeuterated_dists)
  
  compare = merge(theoretical_spectra,
                  observed_spectra,
                  by = c("Peptide", "Time", "IntDiff"))
  
  weights = weights[, colnames(weights) != "Mass", with = FALSE]
  compare = merge(compare, weights, by = c("Peptide", "Time", "IntDiff"))
  compare = compare[ExpectedPeak > 0]
  to_optim_theta = function(theta) {
    compare[, sum((Weight^(2*theta))*(ExpectedPeak - Intensity)^2, na.rm = T)]
  }
  optimized_theta = optim(0.1, to_optim_theta, method = "BFGS")
  optimized_theta$par
}

#' @keywords internal
getExpectedPeaksDerivatives = function(observed_spectra, peptides_cluster,
                                       pept_seg_struct, num_parameters,
                                       time_0_data, undeuterated_dists) {
  time_points = unique(observed_spectra$Time)
  observed_spectra = observed_spectra
  peptides_cluster = peptides_cluster
  pept_seg_struct = pept_seg_struct
  num_parameters = num_parameters
  time_0_data = time_0_data
  undeuterated_dists = undeuterated_dists
  num_ex_peptide = peptides_cluster[, .(Peptide, Segment, MaxUptake)]
  num_ex_peptide = num_ex_peptide[, .(NumEx = sum(MaxUptake)), by = "Peptide"]
  
  function(parameters) {
    analytical_gradient_probs = get_analytical_derivatives(peptides_cluster, pept_seg_struct, undeuterated_dists, num_parameters, time_points, parameters)
    theoretical_spectra = getExpectedSpectra(parameters, pept_seg_struct, num_parameters + 1, observed_spectra, undeuterated_dists)
    compare = merge(observed_spectra, theoretical_spectra,
                    by = c("Peptide", "Time", "IntDiff"),
                    sort = FALSE)
    
    derivatives = lapply(unique(analytical_gradient_probs$Time), function(time) {
      lapply(unique(analytical_gradient_probs$Segment[analytical_gradient_probs$Time == time]), function(segment) {
        lapply(unique(analytical_gradient_probs$ParameterID[analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment]), function(param_id) {
          lapply(unique(analytical_gradient_probs$Peptide[analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment]), function(peptide) {
            undeuterated_probabilities_current = undeuterated_dists[[peptide]]
            num_exchangeable = num_ex_peptide$NumEx[num_ex_peptide$Peptide == peptide]
            
            lapply(observed_spectra[observed_spectra$Peptide == peptide & observed_spectra$Time == time, unique(Charge)], function(charge) {
              total = sum(observed_spectra$Intensity[observed_spectra$Peptide == peptide & observed_spectra$Time == time & observed_spectra$Charge == charge], na.rm = T)
              lapply(1:(length(undeuterated_probabilities_current) + num_exchangeable), 
                     function(k) {
                       peak_range = max(1, k - num_exchangeable):min(length(undeuterated_probabilities_current), k)
                       range_len = length(peak_range)
                       list(Peptide = peptide,
                            Segment = segment,
                            Time = time,
                            Charge = charge,
                            ParameterID = param_id,
                            DerivativeEt = sum(vapply(peak_range,
                                                      function(d) {
                                                        filt = analytical_gradient_probs$ParameterID == param_id & analytical_gradient_probs$Time == time & analytical_gradient_probs$Segment == segment & analytical_gradient_probs$Peptide == peptide & analytical_gradient_probs$NumExchanged == k - d
                                                        x = undeuterated_probabilities_current[d] * total * analytical_gradient_probs$ForDerivative[filt]
                                                        x
                                                      }, numeric(1)), na.rm = TRUE),
                            IntDiff = k - 1)
                     })
            })
          })
        })
      })
    })
    derivatives_df = rbindlist(unlist(unlist(unlist(unlist(unlist(derivatives, F, F), F, F), F, F), F, F), F, F))
    
    res = merge(derivatives_df, compare, 
                by = c("Time", "Peptide", "Charge", "IntDiff"),
                allow.cartesian = T, sort = FALSE)
    res
  }
}

#' Get expected exchange probabilities at the peptide level
#' 
#' @param spectra xxx
#' @param peptides_cluster xxx
#' @param time_0_data xxx
#' @param times xxx
#' @param betas xxx
#' 
#' @export
#' 
getExpectedProbabilitiesTable = function(spectra,
                                         peptides_cluster, 
                                         time_0_data,
                                         times, betas) {
  segments = unique(peptides_cluster[, .(Segment, Start, End, MaxUptake)])
  peptides_cluster[, Present := 1]
  ps_m = data.table::dcast(peptides_cluster, Peptide ~ Segment, value.var = "Present", fill = 0)
  ps_m = ps_m[, c("Peptide", levels(segments$Segment)), with = FALSE]
  num_parameters = segments[["MaxUptake"]]
  
  expected_peak_heights = lapply(
    as.character(unique(peptides_cluster[["Peptide"]])),
    function(peptide) {
      present_segments = as.logical(ps_m[Peptide == peptide, -1, with = FALSE])
      segment_probabilities = makeSegmentProbabilities(betas, times, num_parameters + 1, present_segments)
      
      num_exchangeable = sum(num_parameters[present_segments])
      lapply(seq_along(times), function(i) {
        time = times[i]
        probs = segment_probabilities[[i]]
        
        peptide_probabilities = getExchangeProbabilities(probs, num_exchangeable, FALSE)
        
        lapply(spectra[Peptide == peptide & Time == time, unique(Charge)], function(charge) {
          list(Peptide = peptide, 
               Time = time, 
               Charge = charge,
               NumExchanged = peptide_probabilities$NumExchanged,
               Probability = peptide_probabilities$Probability)
        })
      })
    })
  data.table::rbindlist(unlist(unlist(expected_peak_heights, F, F), F, F))
}

#' Expected exchange probabilities at the segment level
#' 
#' @param spectra xxx
#' @param peptides_cluster xxx
#' @param time_0_data xxx
#' @param times xxx
#' @param betas xxx
#' 
#' @export
#' 
getExpectedSegmentProbabilitiesTable = function(spectra,
                                                peptides_cluster, 
                                                time_0_data,
                                                times, betas) {
  segments = unique(peptides_cluster[, .(Segment, Start, End, MaxUptake)])
  peptides_cluster[, Present := 1]
  ps_m = data.table::dcast(peptides_cluster, Peptide ~ Segment, value.var = "Present", fill = 0)
  ps_m = ps_m[, c("Peptide", levels(segments$Segment)), with = FALSE]
  num_parameters = segments[["MaxUptake"]]
  
  expected_peak_heights = lapply(
    as.character(unique(peptides_cluster[["Segment"]])),
    function(segment) {
      present_segments = colnames(ps_m)[-1] == segment
      segment_probabilities = makeSegmentProbabilities(betas, times, num_parameters + 1, present_segments)
      
      lapply(seq_along(times), function(i) {
        time = times[i]
        probs = segment_probabilities[[i]]
        
        list(Segment = segment, 
             Time = time, 
             NumExchanged = 0:peptides_cluster[Segment == segment, unique(MaxUptake)],
             Probability = probs[[1]])
      })
    })
  data.table::rbindlist(unlist(expected_peak_heights, F, F))
}
