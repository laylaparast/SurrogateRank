# Generate high-dimensional multistudy example data

example.data.highdim.multistudy = generate.example.data.highdim.multistudy(
  epsilon = 0.2,
  M = 5,
  sample_sizes = c(25, 50, 100, 150, 250),
  J = 500,
  prop_valid = 0.1,
  u_tau_min = 0.01,
  u_tau_max = 0.1, 
  u_nu_min = 0.01,
  u_nu_max = 0.1,
  prop_invalid_under = 0.5,
  seed = 12345
)

path = fs::path("data", "example.data.highdim.multistudy.RData")
save(example.data.highdim.multistudy, file = path)

rm(list = ls())
