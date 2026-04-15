# Generate high-dimensional example data

example.data.highdim = generate.example.data.highdim(n1 = 25, 
                                                     n0 = 25, 
                                                     p = 1000,
                                                     prop_valid = 0.1, 
                                                     valid_sigma = 1, 
                                                     corr = 0, 
                                                     mode = "simple",
                                                     seed = 12345
)

path = fs::path("data", "example.data.highdim.RData")
save(example.data.highdim, file = path)

rm(list = ls())

