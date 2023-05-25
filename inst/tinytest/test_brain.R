# Examples verifiable by hand.
# 1. (0.5 + 0.6x) * (0.3 + 0.4x) = 0.15 + x * (0.2 + 0.18) + x ^ 2 * (0.24)

probs_1 = IsoHDX:::getExchangeProbabilities(
  data.table::data.table(Segment = c("A", "B"),
                         Start = c(1, 3),
                         End = c(2, 4),
                         MaxUptake = c(1, 1)),
  by_segment_probabilities = list(c(0.5, 0.6), c(0.3, 0.4)),
  FALSE
)

tinytest::expect_equal(probs_1$Probability, c(0.15, 0.38, 0.24))

# Np. p-stwo segmentu 1, 3 hydrogens: 
# p00 + p01I + p02I + p03I
# równania: