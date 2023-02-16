Social Learning Model

** to switch between density-dependent learning and region-specific learning change for

run_future.py: 
- generates time series for adapters (x), utility (dU), risk (R), mortalities (m) and relative bad air days (BAD)
- requires command line arguments for year, policy, and initial condition
- to run, enter "python run_future.py <year> <policy> <initial condition> <density>"
  or "python3 run_future.py <year> <policy> <initial condition> <density>"
- year can be 2050 or 2100
- policy can be REF, P37 or P45
- initial condition can be IC1, IC2, IC3, IC4 or IC5
- density can be density or no_density
- output is saved in ./output

plot_future.py:
- takes output from run_future.py and plots x, dU, m, and BAD
- to run, enter "python plot_future.py <year> <policy> <initial condition> <density>" 
  or "python3 plot_future.py <year> <policy> <initial condition> <density>"
- density can be density or no_density
- output is saved in ./figures
