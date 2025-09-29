# ✅ Run Simulation Over Parameter Grid
start_time <- Sys.time()

simulation_results <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP))
}))

end_time <- Sys.time()
