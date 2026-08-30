function cfg = config()
% config.m  --  Central configuration for the frontotemporal connectivity pilot
% Dataset: OpenNeuro ds004117 (Sternberg verbal working memory, hd-EEG).
% VALUES BELOW ARE NOW SET from the s00_inventory output of this dataset.
%
% For run order and what each script (s00-s16) produces, see README.md.
% Every script begins with:  cfg = config();

% --------------------- PATHS (EDIT to your machine) ----------------------
cfg.bids_dir   = fullfile('..','data','ds004117');           % EDIT to your local path
cfg.deriv_dir  = fullfile('..','derivatives','ftc');         % EDIT to your local path
cfg.eeglab_dir = fullfile('..','..','eeglab2023.0');         % EDIT to your EEGLAB install

% ---------------- DATASET SUBSET (for a fast first pass) ----------------
% []  = process ALL datasets (full run, may take hours: ICA x ~85 datasets).
% N   = process only the first N datasets to get a quick preliminary figure,
%       then set back to [] for the final run.
cfg.subset_n = [];

% ---------------- PREPROCESSING ----------------
cfg.resample_hz = 250;   % data already 250 Hz; resample is a no-op safeguard
cfg.hp_hz       = 1.0;
cfg.lp_hz       = 45;
cfg.line_hz     = [];    % 50 (EU) if you see a line peak; [] = skip cleanline
cfg.amp_thresh  = 150;   % uV; per-trial maintenance-window rejection

% ---------------- EVENTS  (resolved from s00 inventory) ----------------
% Maintenance/retention onset = the 'show_dash' / work_memory event,
% which appears in the EEG.event.type field as 'WM'.
cfg.ev_maint_lock = 'WM';
cfg.maint_window  = [0.5 2.5];      % s rel. to dash onset; delay is ~3.58 s, so
                                    % this ends ~1 s before the probe. Safe.
% Memory load = number of to-remember letters, stored in event field memory_cond.
cfg.load_field    = 'memory_cond';
cfg.load_low      = {'3'};          % set size 3
cfg.load_high     = {'7'};          % set size 7  (5 = medium, excluded here)

% Trial inclusion from the consistent semantic field task_role (NOT the
% remapped 'type' field, which is not uniform across subjects).
cfg.correct_only            = true;          % keep only correct-response trials
cfg.role_field              = 'task_role';
cfg.role_feedback_correct   = 'feedback_correct';
cfg.role_feedback_incorrect = 'feedback_incorrect';
cfg.role_bad                = 'bad_trial';   % always excluded

% Analyse band power in log space (power is log-normal; a priori, not p-hacking).
cfg.log_power = true;

% ---------------- ROIs (case-insensitive; montage uses UPPERCASE labels) ------
cfg.roi_frontal  = {'F3','Fz','F4','FCz','AFz'};
cfg.roi_temporal = {'T7','T8','TP7','TP8','FT7','FT8'};

% ---------------- BANDS ----------------
cfg.bands.theta = [4 7];
cfg.bands.beta  = [13 30];

% ---------------- CONNECTIVITY ----------------
cfg.filt_order_factor = 3;

% ---------------- STATS ----------------
cfg.stat.n_perm = 5000;
cfg.stat.alpha  = 0.05;
cfg.stat.tail   = 0;

% ---------------- FIGURES ----------------
cfg.fig.dpi = 300;
end
