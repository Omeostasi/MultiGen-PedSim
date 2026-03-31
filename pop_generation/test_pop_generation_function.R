source("pop_generation/pop_generation_function.R")


# TEST 1

founder <- quickHaplo(nInd = 10^6, nChr = 2, segSites = 100)

SP <- SimParam$new(founder)
SP$setSexes("yes_rand")


pop <- pop_generation(founder)

pop@id

pop@mother

pop@father
