# 1 peptide, 1 spectrum, prob(o) = 1: same spectrum is generated
times = 1
ps_m = matrix(c(1),
              byrow = T, nrow = 1, ncol = 1)
ps_m = data.table::as.data.table(ps_m)


segments = "DK"
peptides = "DK"

ps_m = cbind(Peptide = peptides, ps_m)
colnames(ps_m)[2] = segments
num_parameters = 1

betas = runif(sum(num_parameters + 1), -2e-1, -1e-1)
betas2 = betas

undeuterated_dists = lapply(peptides, function(x) BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(x), nrPeaks = 3)$isoDistr)
undeuterated_dists[[1]] = undeuterated_dists[[1]] / sum(undeuterated_dists[[1]])

names(undeuterated_dists) = peptides

seg_probs = IsoHDX:::getSegmentProbabilitiesFromParams(list(c(0, 0)), 1)

observed_spectra_fake = data.table::data.table(Time = 1,
                                               Charge = 2,
                                               Intensity = 1,
                                               Merge = T)

observed_spectrum = IsoHDX:::getExpectedSpectra(c(0, 0), ps_m, 2, observed_spectra = observed_spectra_fake, undeuterated_dists = undeuterated_dists)
observed_spectrum[, IntDiff := 0:3]

# probabilities: 0.5, 0.5
# Intensities according to the model:
# Monoisotopic: nothing was exchanged: 0.5 probability: undeuterated_dists[[1]][1] * 0.5
# 2nd peak: monoisotopic + 1 or 1 + nothing: both 0.5 prob: -> sum(undeuterated_dists[[1]][1:2]*0.5)
# 3rd peak: 2 + 1 exchanged or 3 + no exchanged: sum(undeuterated_dists[[1]][2:3]*0.5)
# 4th peak: 3 + 1 exchanged: undeuterated_dists[[1]][3] * 0.5

tinytest::expect_equal(observed_spectrum$ExpectedPeak,
                       c(undeuterated_dists[[1]][1] * 0.5,
                         sum(undeuterated_dists[[1]][1:2]*0.5),
                         sum(undeuterated_dists[[1]][2:3]*0.5),
                         undeuterated_dists[[1]][3] * 0.5))

# Noiseless optimization

times = c(1, 2)
ps_m = matrix(c(1),
              byrow = T, nrow = 1, ncol = 1)
ps_m = data.table::as.data.table(ps_m)


segments = "DK"
peptides = "DK"

ps_m = cbind(Peptide = peptides, ps_m)
colnames(ps_m)[2] = segments
num_parameters = 1

betas = runif(sum(num_parameters + 1), -2e-1, -1e-1)

undeuterated_dists = lapply(peptides, function(x) BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(x), nrPeaks = 3)$isoDistr)
undeuterated_dists[[1]] = undeuterated_dists[[1]] / sum(undeuterated_dists[[1]])

names(undeuterated_dists) = peptides

seg_probs = IsoHDX:::getSegmentProbabilitiesFromParams(list(c(-0.5, -0.2)), c(1, 2))

observed_spectra_fake = data.table::data.table(Time = 1:2,
                                               Charge = 2,
                                               Intensity = 1,
                                               Merge = T)

observed_spectrum = IsoHDX:::getExpectedSpectra(c(-0.5, -0.2), ps_m, 2, observed_spectra = observed_spectra_fake, undeuterated_dists = undeuterated_dists)
observed_spectrum[, IntDiff := rep(0:3, 2)]

peptide_segment_structure = ps_m
peptides_cluster = data.table::melt(peptide_segment_structure,
                                    id.vars = "Peptide", value.name = "Present", variable.name = "Segment")[Present == 1][, .(Peptide, Segment, Start = NA, End = NA, MaxUptake = 1, Present)]
peptides_cluster$egment = factor(peptides_cluster$Segment, levels = segments, ordered = T)
peptides_cluster$Peptide = factor(peptides_cluster$Peptide, levels = peptides, ordered = T)

time_0_data = data.table::data.table(Peptide = "DK", 
                                     Charge = 2,
                                     NumPeaks = 3)

num_segments = ncol(ps_m) - 1
data.table::setnames(observed_spectrum, "ExpectedPeak", "Intensity")
observed_spectrum$Charge = 2

betas_init_sim_1 = runif(sum(num_parameters) + num_segments, -1e-1, -1e-2)
res_selected_start_sim_1 = IsoHDX::fitIsoSegmentModel(observed_spectrum, peptides_cluster, time_0_data, 
                                                      undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)

tinytest::expect_equal(res_selected_start_sim_1$OptimizationResult$par,
                       c(-0.5, -0.2), tolerance = 1e-4)

observed_spectrum2 = data.table::copy(observed_spectrum)
observed_spectrum2$Intensity = observed_spectrum2$Intensity * (1 + 0.05 * rnorm(8, 0, 1))

tinytest::expect_silent({
  res_selected_start_sim_2 = IsoHDX::fitIsoSegmentModel(observed_spectrum2, peptides_cluster, time_0_data, 
                                                        undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)
})

# TODO: test these final weights
# res_selected_start_sim_2$OptimizationResult$weights.Weight

# Behavior with a missing peak

observed_spectrum3 = data.table::copy(observed_spectrum)
observed_spectrum3 = observed_spectrum3[-2, ]
tinytest::expect_silent({
  res_selected_start_sim_mispeak = IsoHDX::fitIsoSegmentModel(observed_spectrum3, peptides_cluster, time_0_data, 
                                                              undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)
})

betas_init_sim_1 = runif(sum(num_parameters) + num_segments, -1e-1, -1e-2)
res_selected_start_sim_mispeak = IsoHDX:::fitIsoSegmentModel(observed_spectrum3, peptides_cluster, time_0_data, 
                                                    undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)
res_selected_start_sim_mispeak$OptimizationResult$par

tinytest::expect_equal(unlist(res_selected_start_sim_mispeak$FinalComparison[, .(Diff = mean((Intensity - ExpectedPeak)^2))], F, F),
                       0, tolerance = 1e-4)

ols_fun = IsoHDX:::getOptimizationProblem(observed_spectrum3, 
                                          peptides_cluster, time_0_data, undeuterated_dists, weights = NULL, 
                                          theta = 1)
value_at_true = ols_fun(c(-0.5, -0.2))

comp = merge(observed_spectrum,
             observed_spectrum3,
             by = c("Peptide", "Time", "IntDiff", "Charge"),
             all.x = T)
# comp[, Intensity.y := ifelse(is.na(Intensity.y), 0, Intensity.y)]

value_expected = sum((comp$Intensity.x - comp$Intensity.y)^2, na.rm = TRUE)

tinytest::expect_equal(abs(value_at_true - value_expected),
                       0)
# TODO: test merging with multiple time points, peptides and segments. That will prove the loss function is calculated correctly

# Missing peak at the end of the distribution
observed_spectrum4 = data.table::copy(observed_spectrum)
observed_spectrum4 = observed_spectrum4[-4, ]
tinytest::expect_silent({
  res_selected_start_sim_mispeak = IsoHDX::fitIsoSegmentModel(observed_spectrum3, peptides_cluster, time_0_data, 
                                                              undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)
})

betas_init_sim_1 = runif(sum(num_parameters) + num_segments, -1e-1, -1e-2)
res_selected_start_sim_mispeak_end = IsoHDX:::fitIsoSegmentModel(observed_spectrum4, peptides_cluster, time_0_data, 
                                                                 undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)

res_selected_start_sim_mispeak_end$FinalComparison

ols_fun = IsoHDX:::getOptimizationProblem(observed_spectrum4, 
                                          peptides_cluster, time_0_data, undeuterated_dists, weights = NULL, 
                                          theta = 1)
value_at_true = ols_fun(c(-0.5, -0.2))

comp = merge(observed_spectrum,
             observed_spectrum4,
             by = c("Peptide", "Time", "IntDiff", "Charge"),
             all.x = T)
# comp[, Intensity.y := ifelse(is.na(Intensity.y), 0, Intensity.y)]

value_expected = sum((comp$Intensity.x - comp$Intensity.y)^2, na.rm = T)

tinytest::expect_equal(abs(value_at_true - value_expected), 0)

# 
times = c(1, 5, 10)
ps_m = matrix(c(1, 0, 0,
                1, 1, 0,
                0, 1, 0,
                0, 1, 1,
                0, 0, 1),
              byrow = T, nrow = 5, ncol = 3)
ps_m = t(ps_m)
ps_m = data.table::as.data.table(ps_m)

segments = lapply(1:5, function(x) sample(setdiff(Peptides::aaList(), "P"), 3, replace = T))
segments = sapply(segments, function(x) paste(x, collapse = "", sep = ""))
peptides = sapply(1:nrow(ps_m), function(i) {
  select_inds = as.logical(ps_m[i, ])
  peptide = paste(segments[select_inds], sep = "", collapse = "")
  peptide
})

ps_m = cbind(Peptide = peptides, ps_m)
colnames(ps_m)[2:6] = segments
num_parameters = rep(3, 5)

betas = runif(sum(num_parameters + 1), -2e-1, -1e-1)

undeuterated_dists = lapply(peptides, function(x) BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(x), nrPeaks = 3)$isoDistr)
undeuterated_dists = lapply(undeuterated_dists, function(x) x / sum(x))
names(undeuterated_dists) = peptides
observed_spectra2 = data.table::data.table(Time = times, 
                                           Charge = 2,
                                           Intensity = 1,
                                           Merge = TRUE)
observed_spectra2 = merge(observed_spectra2, data.table::data.table(Peptide = peptides, Merge = TRUE), by = "Merge", allow.cartesian = T)
observed_spectra2 = merge(observed_spectra2, data.table::data.table(IntDiff = 0:8, Merge = TRUE), by = "Merge", allow.cartesian = T)
observed_spectra2[, Intensity := Intensity / sum(Intensity), by = c("Time", "Charge", "Peptide")]

sim_spectra = IsoHDX:::getExpectedSpectra(betas,
                                          ps_m,
                                          num_parameters + 1,
                                          observed_spectra2,
                                          undeuterated_dists)
# sim_spectra[, ExpectedPeak := ExpectedPeak + 0.01 * rnorm(.N, 0, 1), by = c("Time", "Peptide")]
sim_spectra[, Charge := 2]
data.table::setnames(sim_spectra, "ExpectedPeak", "Intensity")

# ggplot(sim_spectra, aes(x = IntDiff, ymin = 0, ymax = Intensity)) +
#   geom_linerange() +
#   facet_grid(Time~Peptide) +
#   theme_bw()

peptide_segment_structure = ps_m
peptides_cluster = data.table::melt(peptide_segment_structure,
                                    id.vars = "Peptide", value.name = "Present", variable.name = "Segment")[Present == 1][, .(Peptide, Segment, Start = NA, End = NA, MaxUptake = 3, Present)]
peptides_cluster[, Segment := factor(Segment, levels = segments, ordered = T)]
peptides_cluster[, Peptide := factor(Peptide, levels = peptides, ordered = T)]

time_0_data = data.table::data.table(Sequence = peptides,
                                     Charge = 2,
                                     NumPeaks = 3)
und_d = undeuterated_dists

sim_spectra2 = data.table::copy(sim_spectra)
num_segments = ncol(ps_m) - 1
betas_init_sim_1 = runif(sum(num_parameters) + num_segments, -1e-1, -1e-2)
res_selected_start_sim_1 = IsoHDX:::fitIsoSegmentModel(sim_spectra2, peptides_cluster, time_0_data, 
                                                       undeuterated_dists, betas_init_sim_1, "PLGLS", FALSE, max_iter = 10)

tinytest::expect_true(unlist(res_selected_start_sim_1$FinalComparison[, .(Diff = mean((ExpectedPeak - Intensity)^2))], F, F) < 1e-3)

tinytest::expect_true(mean((res_selected_start_sim_1$OptimizationResult$par - betas)^2) < 1e-2)

sim_spectra3 = data.table::copy(sim_spectra)
sim_spectra3[, Intensity := Intensity + 1e-3 * Intensity * rnorm(.N, 0, 1), by = c("Time", "Peptide")]
sim_spectra3[, Intensity := Intensity / sum(Intensity), by = c("Time", "Peptide")]

random_start_noise = runif(sum(num_parameters) + num_segments, -1e-1, -1e-2)
close_start = betas + runif(length(betas), -1e-5, 1e-5)

res_rand_noise = IsoHDX:::fitIsoSegmentModel(sim_spectra3, peptides_cluster, time_0_data, 
                                             undeuterated_dists, random_start_noise, "PLGLS", FALSE, max_iter = 10)
res_close_noise = IsoHDX:::fitIsoSegmentModel(sim_spectra3, peptides_cluster, time_0_data, 
                                             undeuterated_dists, close_start, "PLGLS", FALSE, max_iter = 10)

tinytest::expect_true(mean((res_rand_noise$OptimizationResult$par - betas)^2) > mean((res_close_noise$OptimizationResult$par - betas)^2))

