source("pop_generation/pop_generation_function.R")


# TEST 1

founder <- quickHaplo(nInd = 20, nChr = 2, segSites = 100)

SP <- SimParam$new(founder)
SP$setSexes("yes_sys")


pop <- pop_generation(founder, p_new_partner = 0.50, overlapping_fraction = 0, nGenerations_pop = 4)

pop@id

pop@mother

pop@father

df <- data.frame(
  ids = pop@id,
  mothers = pop@mother,
  fathers = pop@father
  
)

df
