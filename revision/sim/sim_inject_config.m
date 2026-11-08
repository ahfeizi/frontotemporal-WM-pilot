function cfg = sim_inject_config()
% sim_inject_config.m -- constants for the Dataset 2 injection power curve
% (revision item R3, second half).
%
% A forward-projected frontal-temporal source pair is injected into the REAL
% preprocessed high-load trials at five coupling strengths, and the battery is
% asked how often it notices. This is a power analysis within high load, NOT a
% re-run of the high-low contrast.
%
% Both arms receive the SAME two sources at the SAME amplitude into the SAME
% real trials; only the coupling differs, and the mix is variance preserving,
% so added power is identical in both arms at every coupling level. The kappa = 0
% level is therefore the false-positive row: the arms then differ only by two
% independent noise realisations of source 2.
%
% Shared machinery lives in P2's revision/sim folder (the forward model, the
% cached CSD operator, the batched inference bridge); nothing is duplicated.

% ======================= PATHS (EDIT for a different machine) ============
% P1's own derivatives: manifest.mat, <key>_clean.set, <key>_trials.mat
cfg.deriv_dir     = 'D:\DataSets\Sternberg Working Memory Dataset\ds004117\derivatives\ftc';   % EDIT
cfg.eeglab_dir    = 'C:\Users\Amirhossein\Documents\eeglab2023.0';                             % EDIT
cfg.fieldtrip_dir = 'C:\Users\Amirhossein\Documents\fieldtrip';                                % EDIT
cfg.stats_dir     = 'C:\Users\Amirhossein\Desktop\Prof. Mulert\files\unified_stats';           % EDIT
% P2's revision/sim, which holds sim_util_path.m, sim_util_csd.m,
% sim_util_stats.m and sim_stats_runner.py.
cfg.shared_sim_dir = ['C:\Users\Amirhossein\Desktop\Prof. Mulert\GitHub\pilot2_final\' ...
                      'frontoparietal-gamma-PLV-dwPLI-pilot\revision\sim'];                    % EDIT

cfg.repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.sim_dir   = fileparts(mfilename('fullpath'));
cfg.out_dir   = fullfile(cfg.sim_dir, 'results');
cfg.cache_dir = fullfile(cfg.sim_dir, 'cache');

% ======================= FORWARD MODEL ==================================
cfg.fwd.headmodel = fullfile(cfg.fieldtrip_dir,'template','headmodel','standard_bem.mat');
cfg.fwd.grid_resolution = 10;      % mm
cfg.fwd.grid_unit       = 'mm';
% The montage is taken from each recording's OWN chanlocs, converted by
% FieldTrip's eeglab2fieldtrip, so no coordinate convention is assumed and the
% 68-versus-69-channel question that the supplied BIDS table raises does not
% arise here at all.

% ======================= DATASET (P1 config.m values) ===================
cfg.srate             = 250;                 % config.m cfg.resample_hz
cfg.maint_window      = [0.5 2.5];           % config.m cfg.maint_window
cfg.filt_order_factor = 3;                   % config.m cfg.filt_order_factor
cfg.bands.theta = [4 7];                     % config.m cfg.battery.bands
cfg.bands.beta  = [13 30];
cfg.wavelet_freqs.theta = 4:1:7;
cfg.wavelet_freqs.beta  = 13:2:29;
cfg.n_cycles = 5;                            % config.m cfg.battery.n_cycles
cfg.roi_frontal  = {'F3','Fz','F4','FCz','AFz'};
cfg.roi_temporal = {'T7','T8','TP7','TP8','FT7','FT8'};
cfg.load_condition = 'high';                 % high load only, as the design states
cfg.min_trials = 5;                          % the battery's >= 5 segments rule

cfg.bands_used = {'theta','beta'};           % beta carries no injection and is
                                              % therefore the specificity check
cfg.target_band = 'theta';                   % the band the coupling is injected in
cfg.measures = {'PLV','dwPLI','ImC'};
cfg.refs     = {'average','csd'};

% ======================= INJECTION ======================================
cfg.site_frontal  = 'F3';                    % design: frontal end of the pair
cfg.site_temporal = 'T7';                    % design: temporal end
cfg.inward_frac   = 0.75;                    % radial fraction at which to take
                                              % the nearest grid point
% Five coupling strengths, evenly spaced up to the scenario simulation's 0.5 so
% the two halves of R3 sit on one scale. kappa = 0 is the false-positive row.
cfg.kappas = [0 0.125 0.25 0.375 0.5];
cfg.phase_offset_deg = 90;                    % lagged pair, away from 0 and 180 deg
% Injected source amplitude as a fraction of each subject's OWN observed ROI
% amplitude in the target band, measured from that subject's real data, so the
% injection is scaled per subject rather than to a global constant.
cfg.inject_frac = 0.5;

% ======================= REPETITIONS / SEEDS ============================
% Halved from 200/40 on 2026-11-06, after sim_inject_probe.m measured this stage
% at 78.4 h simulation + 6.3 h inference. What this stage establishes is the
% SHAPE of the power curve: at what coupling strength each measure begins to
% detect, and whether there is a regime where PLV detects while dwPLI and ImC do
% not. The standard error on a detection rate is 0.035 at 200 repetitions and
% 0.050 at 100. If the scenario-stage result holds, the comparison at a given
% kappa is between roughly 5% and roughly 80% -- a gap of 75 points, which 30
% repetitions would resolve. The precision of the threshold itself is set by the
% spacing of the five kappa levels, not by the repetition count. 38 extra hours
% for 0.015 of standard error is not a trade worth making.
%
% The one place precision genuinely matters is the kappa = 0 row, which is a
% false-positive rate. That evidence already exists and is stronger: P2's S0 ran
% 200 studies per dataset and returned 2-6% in every cell. The kappa = 0 row
% here is an independent confirmation on real recordings rather than the primary
% result, and 100 repetitions still bound it usefully.
cfg.n_rep_hilbert = 100;    % repetitions per coupling level, Hilbert path
cfg.n_rep_wavelet = 20;     % the first 20 also run on the wavelet path
cfg.seed_base = 20261030;
cfg.alpha      = 0.05;      % detection: sign-flip p < alpha
cfg.bf_thresh  = 3;
cfg.tost_bound = 0.4;       % matches util_contrast_stats.TOST_BOUND
cfg.use_parfor = true;

cfg.probe_reps = 2;         % sim_inject_probe.m
cfg.probe_subjects = 3;

% ======================= OBSERVED VALUES (for the summary) ==============
% From the R2 results; printed beside the simulated curve, never an input.
cfg.obs.theta_plv_dz   = 0.326;
cfg.obs.theta_plv_bf01 = 2.17;
cfg.obs.beta_plv_dz_range = [-0.028 0.139];
end
