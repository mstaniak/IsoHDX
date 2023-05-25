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
                                  times, max_peaks = 10,
                                  weights = 1) {
  monoisotopic_masses = getMonoisotopicMasses(peptides_cluster, max_peaks)
  observed_spectra_peakid = getPeakIds(observed_spectra, monoisotopic_masses)
  # previous_peaks <<- 1
  
  function(betas) {
    theoretical_spectra = getExpectedSpectra(observed_spectra,
                                             peptides_cluster,
                                             times, betas, max_peaks)
    theoretical_spectra = data.table::rbindlist(unlist(theoretical_spectra, FALSE, FALSE))
    theoretical_spectra = theoretical_spectra[order(Peptide, Time, Mass)]
    theoretical_spectra[, PeakId := 1:.N,  by = c("Peptide", "Time")]
    data.table::setnames(theoretical_spectra, "Peptide", "Peptide")
    compare = merge(theoretical_spectra, observed_spectra_peakid,
                    by = c("Peptide", "Time", "PeakId"),
                    all.x = T, all.y = T)
    # compare = compare[!is.na(Intensity)]
    # crit = (previous_peaks * (compare[["Intensity"]] - compare[["ExpectedPeak"]])) ^ 2
    # crit = ifelse(is.nan(crit), 0, crit)
    # res = sum(crit, na.rm = TRUE)
    # weights <- ifelse(is.finite(1 / sqrt(compare[["ExpectedPeak"]])),
    #                    1 / sqrt(compare[["ExpectedPeak"]]),
    #                    0)
    # weights <- ifelse(is.nan(weights), 0, weights)
    # previous_peaks <<- weights
    # res
    sum( (weights * (compare[["Intensity"]] - compare[["ExpectedPeak"]])) ^ 2,
         na.rm = TRUE)
  }
}


getPeakIds = function(observed_spectra, monoisotopic_masses) {
  data.table::rbindlist(
    lapply(
      split(observed_spectra, observed_spectra[, .(Peptide, Time)]),
      function(x) {
        mono_info = monoisotopic_masses[Peptide == unique(x$Peptide)]
        num_peaks = mono_info$NumExchangeable + mono_info$UndeuteratedLength
        which_peak = guess_peak(mono_info$Monoisotopic, num_peaks, x$Mass)
        peak_id = which_peak:(nrow(x) + which_peak - 1)
        x = x[order(Mass)]
        x$PeakId = peak_id
        x = merge(x, data.table::data.table(Peptide = mono_info$Peptide,
                                            Time = unique(x$Time),
                                            PeakId = 1:num_peaks),
                  by = c("PeakId", "Peptide", "Time"), all.x = TRUE, all.y = TRUE)
        x[, Mass := ifelse(is.na(Mass), mono_info$Monoisotopic + PeakId - 1, Mass)]
        x[, Peptide := Peptide]
        x
      }
    )
  )
}

guess_peak = function(mono, num_peaks, masses) {
  which.min(abs(min(masses) - (unique(mono) + 0:unique(num_peaks))))
}

getMonoisotopicMasses = function(peptides_cluster, max_peaks = 10) {
  peptides_info = peptides_cluster[, .(NumExchangeable = sum(MaxUptake)), 
                                   by = "Peptide"]
  data.table::rbindlist(
    lapply(
      split(peptides_info, peptides_info$Peptide),
      function(x) {
        iso_dist = BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(x$Peptide), 
                                   nrPeaks = max_peaks)
        iso_dist$isoDistr = iso_dist$isoDistr / sum(iso_dist$isoDistr)
        cbind(x,
              Monoisotopic = iso_dist$monoisotopic,
              UndeuteratedLength = length(iso_dist$isoDistr))
      }))
}


#' Calculate expected isotopic peaks based on fixed values of beta parameters
#' @param spectra data.table of observed spectra
#' @inheritParams getOptimizationProblem
#' @param beta numeric vector of fixed values of beta parameters
#' @keywords internal
getExpectedSpectra = function(spectra,
                              peptides_cluster, 
                              times, betas, max_peaks = 10) {
  segments = unique(peptides_cluster[, .(Segment, Start, End, MaxUptake)])
  peptides_cluster[, Present := 1]
  ps_m = data.table::dcast(peptides_cluster, Peptide ~ Segment, value.var = "Present", fill = 0)
  ps_m = ps_m[, c("Peptide", levels(segments$Segment)), with = FALSE]
  num_parameters = segments[["MaxUptake"]]

  expected_peak_heights = lapply(
    unique(peptides_cluster[["Peptide"]]),
    function(peptide) {
      undeuterated_dist_raw = BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(peptide), nrPeaks = max_peaks)
      undeuterated_dist = undeuterated_dist_raw[["isoDistr"]]
      undeuterated_dist = undeuterated_dist / sum(undeuterated_dist)
      monoisotopic = undeuterated_dist_raw[["monoisotopicMass"]]
      present_segments = as.logical(ps_m[Peptide == peptide, -1, with = FALSE])
      segment_probabilities = IsoHDX:::makeSegmentProbabilities(betas, times, num_parameters, present_segments)
      num_exchangeable = sum(num_parameters[present_segments])
      lapply(seq_along(times), function(i) {
        time = times[i]
        probs = segment_probabilities[[i]]
        peptide_probabilities = IsoHDX:::getExchangeProbabilities(probs, num_exchangeable, FALSE)
        spectrum = spectra[Peptide == peptide & Time == time]
        peak_heights = IsoHDX:::getExpectedPeakHeights(spectrum, 
                                                       peptide_probabilities[["Probability"]], undeuterated_dist, 
                                                       num_exchangeable)
        data.table::data.table(Peptide = peptide, Time = time, 
                               UndeuteratedDist = c(undeuterated_dist, rep(NA, length(peak_heights) - length(undeuterated_dist))), 
                               Mass = seq(monoisotopic, 
                                          monoisotopic + length(peak_heights) - 1, 
                                          by = 1), ExpectedPeak = peak_heights)
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
getExchangeProbabilities = function(by_segment_probabilities,
                                    num_exchangeable,
                                    approximate_root = FALSE) {
  no_exchange_probability = prod(sapply(by_segment_probabilities, function(x) x[1]))
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
        # i = j + 1 <=> j = i - 1
        # q_j = -1/j * sum^{l = 1}^{j}q_{j - l}psi(l) <=> psi(j) = ((-j * q_j) - sum_{l = 1}^{j - 1}q_{j - l}psi(l))/q_0 
        # power_sums[i] = -(i - 1) * coefs[i] / coefs[i - 1]  
        # power_sums[i] = (-(i - 1) * coefs[i] - sum(rev(coefs[2:(i)]) * power_sums[2:(i)])) / coefs[1]
        # power_sums[i] = -(i - 1) * sum(rev(coefs[2:(i - 1)]) * power_sums[2:(i - 1)]) / coefs[1]
        # i = 3 <=> j = 2
        # q_2 = (-1/2) * (q_1 * psi(1) + q_0 * psi(2))
        # q_2 * (-2) = q_1 * psi(1) + q_0 * psi(2)
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
  data.table::data.table(NumExchanged = 0:num_exchangeable,
                         Probability = unlist(probabilities, FALSE, FALSE))
}


#' Get peptide-level probabilities from segment-level probabilities
#' @keywords internal
getExchangeProbabilitiesOld = function(by_segment_probabilities,
                                    num_exchangeable,
                                    approximate_root = FALSE) {
  roots_by_segment = lapply(by_segment_probabilities, function(x) pracma::roots(rev(x)))
  all_roots = unlist(roots_by_segment, FALSE, FALSE)
  no_exchange_probability = prod(sapply(by_segment_probabilities, function(x) x[1]))

  polynomials = vector("list", num_exchangeable)
  roots_list = lapply(seq_along(polynomials), function(x) sum(all_roots ^ (-x)))
  roots_list = sapply(roots_list, function(x) Re(x))
  for (i in seq_len(num_exchangeable)) {
    if (i == 1) {
      polynomials[[i]] = Re(- no_exchange_probability * sum(as.complex(1)/ all_roots))
    } else {
      polynomials[[i]] = Re((-1 / i) * sum(Re(rev(c(no_exchange_probability, unlist(polynomials[1:(i - 1)], FALSE, FALSE)))) * roots_list[1:i]))
    }
  }
  probabilities = lapply(polynomials, Re)
  probabilities = c(no_exchange_probability, probabilities)
  data.table::data.table(NumExchanged = 0:num_exchangeable,
                         Probability = unlist(probabilities, FALSE, FALSE))
}

#' Convert beta parameters to probabilities
#' @keywords internal
getProbabilitiesFromBetas = function(betas, time) {
  total = 1 + sum(exp(betas * time))
  c(exp(betas * time), 1) / total
}
