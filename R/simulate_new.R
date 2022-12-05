#' Simulate spectra from a model that relates links peptide-level spectra to segment-level exchange probabilities
#' @export
simulateSegmentModel = function(peptides_cluster,
                                time_points,
                                beta_parameters,
                                error_sd = 0.1, total = 1e5,
                                max_iso_peaks = 10) {
  undeuterated_dists_raw = getUndeuteratedDists(peptides_cluster,
                                                max_iso_peaks)
  
  getSpectraByTime(peptides_cluster,
                   undeuterated_dists_raw,
                   time_points,
                   beta_parameters,
                   error_sd, total)
}


#' Create simulated cluster of peptides that share segments
#' @export
createPeptidesCluster = function(segment_peptide_matrix, num_exchangeable) {
  if(is.null(colnames(segment_peptide_matrix))) {
    num_cols = ncol(segment_peptide_matrix)
    colnames(segment_peptide_matrix) = paste0("Pept", seq_len(num_cols))
  }
  if (is.null(row.names(segment_peptide_matrix))) {
    segments = paste0("Seg", seq_len(nrow(segment_peptide_matrix)))
  } else {
    segments = row.names(segment_peptide_matrix)
  }
  aa_vec= getAminoAcids()
  aa_vec = setdiff(aa_vec, "P")
  
  segment_peptide = data.table::as.data.table(segment_peptide_matrix)
  segment_peptide[, MaxUptake := num_exchangeable]
  sequences = sapply(seq_len(nrow(segment_peptide)),
                     function(i) paste0(sample(aa_vec, segment_peptide[, MaxUptake][i] + 1),
                                        collapse = ""))
  segment_peptide[, Segment := sequences] #paste0(sequences, 1:.N) 
  pept_cols = setdiff(colnames(segment_peptide), c("Segment", "MaxUptake"))
  segment_peptide = data.table::melt(segment_peptide, id.vars = c("Segment", "MaxUptake"))
  segment_peptide = segment_peptide[value == 1]
  segment_peptide[, value := NULL]
  segment_peptide[, Peptide := paste0(unique(Segment), collapse = ""), by = "variable"]
  segment_peptide[, variable := NULL]
  segment_peptide[, Start := NA_real_]
  segment_peptide[, End := NA_real_]
  segment_peptide
}

getUndeuteratedDists = function(peptides_cluster, max_iso_peaks) {
  peptides = unique(peptides_cluster[, Peptide])
  dists_raw = lapply(peptides, function(peptide) {
    BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(peptide), nrPeaks = max_iso_peaks)
  })
  names(dists_raw) = peptides
  dists_raw = lapply(dists_raw, function(x) {
    x$isoDistr = x$isoDistr / sum(x$isoDistr)
    x
  })
  dists_raw
}

getSpectraByTime = function(peptides_cluster,
                            undeuterated_dists_raw,
                            time_points,
                            beta_parameters,
                            error_sd, total) {
  data.table::rbindlist(lapply(time_points, function(time) {
    getSpectraByPeptide(peptides_cluster,
                        undeuterated_dists_raw,
                        beta_parameters,
                        time,
                        error_sd, total
                        
    )
  }))
}

getSpectraByPeptide = function(peptides_cluster,
                               undeuterated_dists_raw,
                               beta_parameters,
                               time,
                               error_sd, total) {
  peptides_cluster[, Present := 1]
  segment_peptide = data.table::dcast(peptides_cluster, Peptide ~ Segment, 
                                      value.var = "Present", fill = 0)
  segment_peptide = segment_peptide[, c("Peptide", unique(peptides_cluster$Segment)),
                                    with = FALSE]
  nums_exchangeable = peptides_cluster[, .(num_ex = sum(MaxUptake)),
                                       by = "Peptide"]
  
  data.table::rbindlist(
    lapply(unique(peptides_cluster$Peptide),
           function(peptide) {
             undeuterated_dist_raw = undeuterated_dists_raw[[peptide]]
             num_exchangeable = nums_exchangeable[Peptide == peptide, num_ex]
             membership = segment_peptide[Peptide == peptide]
             membership = membership[, colnames(membership) != "Peptide",
                                     with = FALSE]
             membership = unlist(membership, FALSE, FALSE) # todo: function out of this
             spectrum = getSpectrum(undeuterated_dist_raw,
                                    membership,
                                    beta_parameters,
                                    time, num_exchangeable, error_sd, total)
             spectrum[["Peptide"]] = peptide
             spectrum
           }))
}

getSpectrum = function(undeuterated_dist_raw,
                       membership,
                       beta_parameters,
                       time, num_exchangeable, error_sd, total) {
  segment_probs = getTheoreticalProbabilities(beta_parameters[as.logical(membership)],
                                              time)
  probabilities = getExchangeProbabilities(segment_probs, num_exchangeable)
  peaks = getSimulatedPeaks(probabilities, undeuterated_dist_raw[["isoDistr"]],
                            num_exchangeable, total, error_sd)
  masses = getSimulatedMasses(undeuterated_dist_raw,
                              num_exchangeable)
  list(Time = time,
       Mass = masses,
       Intensity = peaks)
}

getTheoreticalProbabilities = function(beta_parameters, time) {
  lapply(beta_parameters, getProbabilitiesFromBetas, time = time)
}

getSimulatedPeaks = function(probabilities,
                             undeuterated_dist,
                             num_exchangeable,
                             total, error_sd) {
  probabilities = probabilities[["Probability"]]
  
  peaks = sapply(1:(length(undeuterated_dist) + num_exchangeable), 
                 function(k) {
                   sum(sapply(max(1, k - num_exchangeable):min(length(undeuterated_dist), k),
                              function(d) {
                                probabilities[k - d + 1] * undeuterated_dist[d] * total
                              }), na.rm = TRUE)
                 })
  peaks = peaks + error_sd * peaks * rnorm(length(peaks), mean = 0, sd = 1)
  peaks
}

getSimulatedMasses = function(undeuterated_dist_raw,
                              num_exchangeable) {
  undeuterated_dist_length = length(undeuterated_dist_raw[["isoDistr"]])
  monoisotopic_mass = undeuterated_dist_raw[["monoisotopicMass"]]
  seq(from = monoisotopic_mass,
      to = monoisotopic_mass + undeuterated_dist_length + num_exchangeable - 1,
      by = 1)  
}

getAminoAcids = function() {
  c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
    "M", "F", "P", "S", "T", "W", "Y", "V") #, "O", "U
} 
