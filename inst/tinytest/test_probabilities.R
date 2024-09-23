# Nowy wzór na prawdopodobieństwa: wymaga dodania do wektorów wejściowych wyrazu wolnego

total = 1 + 2 * exp(1)
tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(1, 1), 1),
                       c(exp(1) / total, exp(1) / total, 1 - (2 * exp(1) / (total))))

total2 = 1 + 2 * exp(2)
tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(c(1, 1), 2),
                       c(exp(2) / total2, exp(2) / total2, 1 - (2 * exp(2) / (total2))))


total3 = 1 + sum(exp(rep(0.5, 4) * 3))
true3 = c(exp(rep(0.5, 4) * 3) / total3, 1 - sum(exp(rep(0.5, 4) * 3) / total3)) 

tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas(rep(0.5, 4), 3),
                       true3)

total4 = 1 + sum(exp((1:5 / 10) * 4))
true4 = c(exp((1:5 / 10) * 4) / total4, 1 - sum(exp((1:5 / 10) * 4) / total4)) 

tinytest::expect_equal(IsoHDX:::getProbabilitiesFromBetas((1:5 / 10), 4),
                       true4)
