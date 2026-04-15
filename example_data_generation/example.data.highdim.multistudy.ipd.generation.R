# Generate high-dimensional example data

example.data.highdim.multistudy.ipd = generate.example.data.highdim.multistudy.ipd(
  M = 5,
  n1 = 25,
  n0 = 25,
  p = 1000,
  prop_valid = 0.1,
  valid_sigma = 1,
  corr = 0,
  mode = "simple",
  seed = 12345
)

path = fs::path("data", "example.data.highdim.multistudy.ipd.RData")
save(example.data.highdim.multistudy.ipd, file = path)

rm(list = ls())
