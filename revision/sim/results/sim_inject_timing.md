# Injection power curve: run-time probe (sim_inject_probe.m)
Measured on this machine, 2026-11-06, over 3 subject(s) x 2 rep(s).
- forward model + CSD operator: 0.2 s, once (cached afterwards)
- loading one subject: 1.5 s (2103619 samples, 96 high-load trials over the probe set)
- one arm, Hilbert path: 7.46 s; one arm, Hilbert+wavelet: 21.28 s
- inference: 1.570 s per cell at N = 23
Counts: 5 coupling levels, 200 reps (40 also on the wavelet path), 6 arms per rep.
**Projection: 78.4 h simulation + 6.3 h inference + 0.6 min loading = 84.7 h total.**
14400 inference cells in all; sim_stats_runner.py parallelises them internally.
The simulation part is a plain subject loop and is not parallelised, because each
subject holds its whole recording in memory. Run it overnight, or lower
cfg.n_rep_hilbert / cfg.n_rep_wavelet, which scale the total linearly.
The linearity shortcut noted at the top of sim_inject_power.m would cut the
simulation part by roughly an order of magnitude if you want it.
