#' Create optimization problem to estimate segment-specific exchange rates via OLS
#' 
#' @param observed_spectra data.table with observed isotopic distributions
#' @param peptides_cluster data.table with information about peptides cluster
#' @param times numeric vector of time points
#' 
#' @export
#' 
getOptimizationProblem = function(observed_spectra,
                                  peptides_cluster,
                                  times) {
  to_optimize = function(betas) {
    theoretical_spectra = getExpectedSpectra(observed_spectra,
                                             peptides_cluster,
                                             times, betas)
    theoretical_spectra = rbindlist(unlist(theoretical_spectra, FALSE, FALSE))
    setnames(theoretical_spectra, "Sequence", "Peptide")
    compare = merge(theoretical_spectra, observed_spectra,
                    by = c("Peptide", "Time", "Mass")) # does not work as expected
    sum(((compare[["Intensity"]] - compare[["ExpectedPeak"]])) ^ 2)
  }
  to_optimize
}


#' Calculate expected isotopic peaks based on fixed values of beta parameters
#' @param spectra data.table of observed spectra
#' @inheritParams getOptimizationProblem
#' @param beta numeric vector of fixed values of beta parameters
#' @keywords internal
getExpectedSpectra = function(spectra,
                              peptides_cluster, 
                              times, betas) {
  segments = unique(peptides_cluster[, .(Segment, Start, End, MaxUptake)])
  peptides_cluster[, Present := 1]
  ps_m = dcast(peptides_cluster, Sequence ~ Segment, value.var = "Present", fill = 0)
  num_parameters = segments[["MaxUptake"]]
  num_time_points = length(times)
  
  expected_peak_heights = lapply(
    unique(peptides_cluster[["Sequence"]]),
    function(peptide) {
      undeuterated_dist_raw = useBRAIN(getAtomsFromSeq(peptide), nrPeaks = 10)
      undeuterated_dist = undeuterated_dist_raw[["isoDistr"]]
      undeuterated_dist = undeuterated_dist / sum(undeuterated_dist) # cheating
      monoisotopic = undeuterated_dist_raw[["monoisotopicMass"]]
      
      present_segments = as.logical(ps_m[Sequence == peptide, -1, with = FALSE]) 
      segment_probabilities = makeSegmentProbabilities(betas,
                                                       times, 
                                                       num_parameters,
                                                       present_segments)
      num_exchangeable = sum(num_parameters[present_segments])
      
      lapply(seq_along(times), function(i) {
        time = times[i]
        probs = segment_probabilities[[i]]
        peptide_probabilities = getExchangeProbabilities(peptides_cluster[Sequence == peptide,
                                                                          .(Segment, Start, End, MaxUptake)],
                                                         probs, FALSE)
        spectrum = spectra[Peptide == peptide & Time == time]
        
        peak_heights = getExpectedPeakHeights(spectrum, 
                                              peptide_probabilities[["Probability"]], 
                                              undeuterated_dist, num_exchangeable)
        
        data.table(Sequence = peptide,
                   Time = time,
                   Mass = seq(monoisotopic, monoisotopic + length(peak_heights) - 1, by = 1),
                   ExpectedPeak = peak_heights)
      })
    })
}


#' Get segment-level exchange probabilities from beta parameters
#' @keywords internal
makeSegmentProbabilities = function(betas, times, num_parameters, present_segments) {
  segment_probs = lapply(times, function(time) {
    all_probs = lapply(seq_along(num_parameters), function(i) {
      if (i == 1) {
        subset_betas = betas[1:num_parameters[i]]
      } else {
        subset_betas = betas[(sum(num_parameters[1:(i - 1)]) + 1):(sum(num_parameters[1:i]))]
      }
      getProbabilitiesFromBetas(subset_betas, time)
    })
    all_probs[present_segments]
  })
}

#' Calculate expected peak heights
#' @keywords internal
getExpectedPeakHeights = function(spectrum,
                                  probabilities,
                                  undeuterated_dist,
                                  num_exchangeable) {
  total = spectrum[, sum(Intensity)]
  probabilities = peptide_probabilities$Probability
  spectrum = spectrum[order(Mass)]
  
  sapply(1:(length(undeuterated_dist) + num_exchangeable), 
         function(k) {
           sum(sapply(max(1, k - num_exchangeable):min(length(undeuterated_dist), k),
                      function(d) {
                        probabilities[k - d + 1] * undeuterated_dist[d] * total
                      }), na.rm = TRUE)
         })
}

#' Get peptide-level probabilities from segment-level probabilities
#' @keywords internal
getExchangeProbabilities = function(segment_data,
                                    by_segment_probabilities,
                                    approximate_root) {
  roots_by_segment = lapply(by_segment_probabilities, function(x) pracma::roots(rev(x)))
  all_roots = unlist(roots_by_segment, FALSE, FALSE)
  no_exchange_probability = prod(sapply(by_segment_probabilities, function(x) x[1]))
  num_exchangeable = sum(segment_data[, MaxUptake])
  
  polynomials = vector("list", num_exchangeable)
  roots_list = lapply(seq_along(polynomials), function(x) sum(all_roots ^ (-x)))
  roots_list = sapply(roots_list, function(x) Re(x))
  for (i in seq_len(num_exchangeable)) {
    if (i == 1) {
      polynomials[[i]] = - no_exchange_probability * sum(as.complex(1)/ all_roots)
    } else {
      polynomials[[i]] = (-1 / i) * sum(rev(c(no_exchange_probability, unlist(polynomials[1:(i - 1)], FALSE, FALSE))) * roots_list[1:i])
    }
  }
  probabilities = lapply(polynomials, Re)
  probabilities = c(no_exchange_probability, probabilities)
  data.table(NumExchanged = 0:num_exchangeable,
             Probability = unlist(probabilities, FALSE, FALSE))
}

#' Convert beta parameters to probabilities
#' @keywords internal
getProbabilitiesFromBetas = function(betas, time) {
  total = 1 + sum(exp(betas * time))
  c(exp(betas * time), 1) / total
}
