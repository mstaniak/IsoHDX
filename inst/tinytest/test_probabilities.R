total = 1 + 2 * exp(1)
tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(0, 1, 1), 1),
                       c(exp(1) / total, exp(1) / total, 1 - (2 * exp(1) / (total))))

total2 = 1 + 2 * exp(2)
tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(0, 1, 1), 2),
                       c(exp(2) / total2, exp(2) / total2, 1 - (2 * exp(2) / (total2))))


total3 = 1 + sum(exp(rep(0.5, 4) * 3))
true3 = c(exp(rep(0.5, 4) * 3) / total3, 1 - sum(exp(rep(0.5, 4) * 3)) / total3)

tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(0, rep(0.5, 4)), 3),
                       true3)

total4 = 1 + sum(exp((1:5 / 10) * 4))
true4 = c(exp((1:5 / 10) * 4) / total4, 1 - sum(exp((1:5 / 10) * 4) / total4)) 

tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(0, (1:5 / 10)), 4),
                       true4)

# non-zero intercept
# time: 4
# intercept: 0.5
# 0: exp(0.1 * 4 - 0.5) / total
# 1: exp(0.2 * 4 - 0.5) / total
# 2: 1 - (exp(0.1 * 4 - 0.5) / total + exp(0.2 * 4 - 0.5) / total)
total_nz = 1 + exp(0.1 * 4 + 0.5) + exp(0.2 * 4 + 0.5)
probs_nz = c(exp(0.1 * 4 + 0.5) / total_nz, exp(0.2 * 4 + 0.5) / total_nz)
probs_nz = c(probs_nz, 1 - sum(probs_nz))

tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(0.5, 0.1, 0.2), 4),
                       probs_nz)

