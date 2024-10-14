# Examples verifiable by hand.
# 1. (0.5 + 0.6x) * (0.3 + 0.4x) = 0.15 + x * (0.2 + 0.18) + x ^ 2 * (0.24)

# data.table::data.table(Segment = c("A", "B"),
#                        Start = c(1, 3),
#                        End = c(2, 4),
#                        MaxUptake = c(1, 1)),

probs_1 = IsoHDX:::getExchangeProbabilities(
  by_segment_probabilities = list(c(0.5, 0.6), c(0.3, 0.4)), 2,
  FALSE
)

tinytest::expect_equal(probs_1$Probability, c(0.15, 0.38, 0.24))

# (p10 + p11I + p12I^2)(p20 + p21I + p22I^2)
# = (p10p20) + I(p10p21 + p11p20) + I^2(p10p22 + p20p12 + p11p21) + I^3(p11p22 + p12p21) + I^4(p12p22)
#

# All non-zero
# 0.2 + I(0.16 + 0.15) + I^2(0.04 + 0.1 + 0.12) + I^3(0.03 + 0.08) + I^4(0.02)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.4, 0.3, 0.2), c(0.5, 0.4, 0.1)), 4)$Probability,
                       c(0.2, 0.31, 0.26, 0.11, 0.02))
# Leading zeros
# 0 + I(0.15 + 0) + I^2(0 + 0.35 + 0.12) + I^3(0.03 + 0.28) + I^4(0.07)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.0, 0.3, 0.7), c(0.5, 0.4, 0.1)), 4)$Probability,
                       c(0.0, 0.15, 0.47, 0.31, 0.07))
# Trailing zeros
# 0.3 + I(0.24 + 0.2) + I^2(0.06 + 0 + 0.16) + I^3(0.04 + 0) + I^4(0)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.6, 0.4, 0.0), c(0.5, 0.4, 0.1)), 4)$Probability,
                       c(0.3, 0.44, 0.22, 0.04, 0.0))

# 0.12 + I(0.28 + 0.18) + I^2(0.42) + I^3(0) + I^4(0)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.4, 0.6, 0.0), c(0.3, 0.7, 0.0)), 4)$Probability,
                       c(0.12, 0.46, 0.42, 0.0, 0.0))
# Zeros in the middle
# 0.2 + I(0.16 + 0) + I^2(0.04 + 0 + 0.3) + I^3(0 + 0.24) + I^4(0.06)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.4, 0.0, 0.6), c(0.5, 0.4, 0.1)), 4)$Probability,
                       c(0.2, 0.16, 0.34, 0.24, 0.06))
# Leading and trailing zeros
# 0 + I(0.12 + 0) + I^2(0.18 + 0 + 0.28) + I^3(0.42 + 0) + I^4(0)
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.3, 0.7, 0.0), c(0.0, 0.4, 0.6)), 4)$Probability,
                       c(0.0, 0.12, 0.46, 0.42, 0.0))

# Small values
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(1 - 1e-6, 1e-6), c(1e-8, 1 - 1e-8)), 2)$Probability,
                       c(0, 1, 0), tolerance = 1e-5)
tinytest::expect_equal(sum(IsoHDX:::getExchangeProbabilities(list(c(1 - 1e-6, 1e-6), c(1e-8, 1 - 1e-8)), 2)$Probability), 1)

tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(1 - 1e-8, 1e-8), c(1e-10, 1 - 1e-10)), 2)$Probability,
                       c(0, 1, 0), tolerance = 1e-5)
tinytest::expect_equal(sum(IsoHDX:::getExchangeProbabilities(list(c(1 - 1e-8, 1e-8), c(1e-10, 1 - 1e-10)), 2)$Probability), 1)

# Longer polynomials [=larger segments]
tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(c(0.3, 0.7, 0.0), c(0.0, 0.4, 0.6)), 4)$Probability,
                       c(polynom:::coef.polynomial(polynom::polynomial(c(0.3, 0.7, 0.0)) * polynom::polynomial(c(0.0, 0.4, 0.6))), 0.0))


seg_1 = c(0.0, 0.5, 0.4, 0.1)
seg_2 = c(0.1, 0.3, 0.3, 0.2, 0.1)
seg_3 = c(0.2, 0.4, 0.4, 0.0)
seg_4 = c(0.1, 0.3, 0.4, 0.2, 0.0)
n = 14

IsoHDX:::getExchangeProbabilities(list(seg_1, seg_2, seg_3, seg_4), 14)

p1 = polynom::polynomial(seg_1)
p2 = polynom::polynomial(seg_2)
p3 = polynom::polynomial(seg_3)
p4 = polynom::polynomial(seg_4)

prod = p1 * p2 * p3 * p4


polynom_four_segments = polynom:::coef.polynomial(prod)

tinytest::expect_equal(IsoHDX:::getExchangeProbabilities(list(seg_1, seg_2, seg_3, seg_4), 14)$Probability,
                       c(polynom_four_segments, 0.0, 0.0))

