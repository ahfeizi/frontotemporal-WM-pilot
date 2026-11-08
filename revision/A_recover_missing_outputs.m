% A_recover_missing_outputs.m
% Recovers A_ppc_full_trials.csv and A_trial_equalisation.md
% without re-running the full 47-subject x 200-rep loop.
%
% Priority:
%   1. Load ppc_full from A_workspace.mat (if it exists in revision/).
%   2. Note that tmp_A_contrasts.csv cannot substitute (it holds X_sub,
%      not ppc_full; and it was deleted after Python succeeded anyway).
%   3. If ppc_full is unavailable, report exactly what to do and still
%      write A_trial_equalisation.md (which is static text).

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));   % adds repo root so config() is found

N      = 47;
N_REPS = 200;

% Band names must match the order config() uses (which is what the original
% run used via BANDS = fieldnames(cfg.battery.bands)).
try
    cfg   = config();
    BANDS = fieldnames(cfg.battery.bands);
catch
    % Fallback: order read from A_trial_equalisation.csv row order.
    BANDS = {'theta','beta','beta_low','beta_high'};
    warning('A_recover:no_config', ...
        'config() not found; using BANDS from CSV row order. Verify this matches the original run.');
end

REFS    = {'average','csd'};
MEAS    = {'PLV','dwPLI','ImC'};
DECOMPS = {'hilbert','wavelet'};

% Reconstruct cell_defs in the exact same loop order as A_trial_equalisation.m
% (band > decomp > ref > measure) so column indices of ppc_full match.
cell_defs = {};
for b = 1:numel(BANDS)
    for d = 1:2, for r = 1:2, for m = 1:3
        cell_defs{end+1} = struct('band',BANDS{b},'measure',MEAS{m},...
            'decomposition',DECOMPS{d},'reference',REFS{r}); %#ok<SAGROW>
    end, end, end
end
n_cells = numel(cell_defs);   % expected: 48
fprintf('cell_defs: %d cells (expected 48)\n', n_cells);

ppc_full = [];
source   = '';

% ── Option 1: workspace dump ─────────────────────────────────────────────
ws_file = fullfile(repo, 'A_workspace.mat');
if isfile(ws_file)
    fprintf('Loading ppc_full from A_workspace.mat ...\n');
    W = load(ws_file, 'ppc_full');
    if isfield(W, 'ppc_full')
        ppc_full = W.ppc_full;
        source   = 'A_workspace.mat';
        fprintf('  Loaded ppc_full [%s]\n', num2str(size(ppc_full)));
    else
        fprintf('  WARNING: A_workspace.mat exists but ppc_full variable is absent.\n');
    end
else
    fprintf('A_workspace.mat not found in revision/.\n');
end

% ── Option 2: tmp files ───────────────────────────────────────────────────
if isempty(ppc_full)
    tmp_in = fullfile(repo, 'tmp_A_contrasts.csv');
    if isfile(tmp_in)
        fprintf('tmp_A_contrasts.csv found, but it contains X_sub (subsampled\n');
        fprintf('contrasts per rep), NOT ppc_full (full-trial PPC). Cannot reconstruct\n');
        fprintf('ppc_full from it.\n');
    else
        fprintf('tmp files not found — deleted by A_trial_equalisation.m after Python\n');
        fprintf('succeeded (expected; main CSV is valid).\n');
    end
end

% ── Option 3: report and advise ──────────────────────────────────────────
if isempty(ppc_full)
    fprintf('\n=== ppc_full unavailable: A_ppc_full_trials.csv CANNOT be written. ===\n\n');
    fprintf('Two recovery paths (choose one):\n\n');
    fprintf('  (a) Minimal re-run — add ONE line to A_trial_equalisation.m immediately\n');
    fprintf('      AFTER line 141 (end of the PPC loop, before the subsampling block):\n\n');
    fprintf('          save(fullfile(repo,''A_workspace.mat''),''ppc_full'',''cell_defs'');\n\n');
    fprintf('      Re-run A_trial_equalisation.m (the 200-rep loop will re-run too;\n');
    fprintf('      to avoid that, comment out lines 143-197 before this run).\n');
    fprintf('      Then re-run this script.\n\n');
    fprintf('  (b) Extract-and-run — copy lines 63-141 of A_trial_equalisation.m\n');
    fprintf('      (subject load + PPC loop only, no subsampling) into a new script\n');
    fprintf('      A_ppc_only.m, add the save() call above, and run A_ppc_only.m.\n');
    fprintf('      That takes minutes, not hours.\n\n');
else
    % ── Write A_ppc_full_trials.csv ──────────────────────────────────────
    ppc_rows = {};
    for c = 1:n_cells
        cd_c = cell_defs{c};
        if ~strcmp(cd_c.measure,'PLV') || ~strcmp(cd_c.decomposition,'hilbert'), continue; end
        x    = ppc_full(:, c);
        x    = x(isfinite(x));
        if numel(x) < 3, continue; end
        mn   = mean(x);  sd_x = std(x);  n_pp = numel(x);
        dz_pp = mn / sd_x;
        t_pp  = mn / (sd_x / sqrt(n_pp));
        ppc_rows{end+1} = {cd_c.band, cd_c.reference, n_pp, dz_pp, t_pp}; %#ok<SAGROW>
    end
    if ~isempty(ppc_rows)
        ppc_T = cell2table(vertcat(ppc_rows{:}), 'VariableNames', ...
            {'band','reference','n','dz_ppc_contrast','t_ppc_contrast'});
        writetable(ppc_T, fullfile(repo, 'A_ppc_full_trials.csv'));
        fprintf('Wrote A_ppc_full_trials.csv (%d rows, source: %s)\n', height(ppc_T), source);
    else
        fprintf('No PLV/hilbert cells with n >= 3. A_ppc_full_trials.csv not written.\n');
    end
end

% ── Always write A_trial_equalisation.md (static content) ────────────────
md = {};
md{end+1} = '# A. Trial-count equalisation';
md{end+1} = '';
md{end+1} = sprintf('N = %d subjects; trial counts: low 9–44, high 11–39.', N);
md{end+1} = sprintf('Subsample to min(n_low, n_high) per participant, %d reps (seeds 1–%d).', N_REPS, N_REPS);
md{end+1} = 'Seed per (rep, subject) = rep × 10^6 + subject_index.';
md{end+1} = 'Wavelet: subsampled at trial level per frequency bin, averaged within band.';
md{end+1} = 'PPC: trial-count-unbiased PLV^2 = (n*PLV^2 - 1)/(n-1), full trials, PLV/Hilbert.';
md{end+1} = '';
md{end+1} = 'Outputs:';
md{end+1} = '  A_trial_equalisation.csv: 48 rows; median dz, 2.5/97.5 pct, median BF01.';
md{end+1} = '  A_ppc_full_trials.csv:    PPC x {average, CSD} x 4 bands = 8 rows.';
fid = fopen(fullfile(repo, 'A_trial_equalisation.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Wrote A_trial_equalisation.md\n');
