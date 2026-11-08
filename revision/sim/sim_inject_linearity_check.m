% sim_inject_linearity_check.m -- Verify the linearity shortcut in sim_inject_power.m.
%
% Picks the first available subject and first run, generates one set of injection
% sources, then compares TWO pipelines:
%
%   FULL: avg_ref(real + inj) -> CSD(real + inj) -> ROI rows
%   SPLIT: [avg_ref(real) + avg_ref(inj)] -> [CSD(real) + CSD(inj)] -> ROI rows
%
% Both must produce identical ROI rows (max abs difference < 1e-9 * signal range).
%
% Run once before sim_inject_power.  Must print PASSED for all checks.
% Expected run time: <30 s (one subject, one run, no repetition loop).
%
% Run:  cd <P1>/revision/sim; sim_inject_linearity_check

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = sim_inject_config();
repo = cfg.repo_root;
addpath(repo);
addpath(cfg.shared_sim_dir);
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end
if ~isfolder(cfg.cache_dir), mkdir(cfg.cache_dir); end

sim_util_path(cfg, {'ft_defaults','ft_prepare_leadfield','ft_prepare_sourcemodel', ...
                    'ft_scalpcurrentdensity','ft_convert_units','ft_datatype_sens'});

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
subj = char(subjects(1));
fprintf('Linearity check subject: %s\n', subj);

keys = manifest.key(manifest.subject == subjects(1));
key = char(keys(1));
EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;

have = {EEG.chanlocs.labels};
labF = cfg.roi_frontal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_frontal));
labT = cfg.roi_temporal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_temporal));
roiLabs = [labF labT];

ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
F = sim_inject_forward(cfg, ft_raw.elec, ft_raw.label, 'ds004117');
[Cfull, csd_labels] = sim_util_csd(cfg, F.elec, F.labels, cfg.srate, 'ds004117');
iCsd = local_idx(csd_labels, roiLabs, 'CSD output');
Crow = Cfull(iCsd, :);
iFwd = local_idx(F.labels, roiLabs, 'forward montage');
gF = F.sites.frontal.gain;
gT = F.sites.temporal.gain;

band = cfg.bands.(cfg.target_band);
fc   = mean(band);
tau  = (cfg.phase_offset_deg/360) * cfg.srate / fc;

rng(42, 'twister');
X_real = double(EEG.data);
n = size(X_real, 2);
sA = util_bandpass(randn(1,n), cfg.srate, band, cfg.filt_order_factor);
sA = sA / std(sA);
sA_del = local_delay(sA, tau);
kap = 0.5;
s2 = kap * sA_del + sqrt(max(1-kap^2,0)) * util_bandpass(randn(1,n), cfg.srate, band, cfg.filt_order_factor);
s2 = s2 / std(s2);

obs_roi_sd = mean(std(X_real(iFwd,:), 0, 2));
sF_sc = cfg.inject_frac * obs_roi_sd / mean(abs(gF(iFwd)));
sT_sc = cfg.inject_frac * obs_roi_sd / mean(abs(gT(iFwd)));

inj_full = (sF_sc*gF) * sA + (sT_sc*gT) * s2;   % nCh x T

% ---- FULL PIPELINE: what the original code did, and the ground truth ----
X_comb   = X_real + inj_full;
X_ar     = X_comb - mean(X_comb, 1);             % average reference of combined
full_avg  = X_ar(iFwd, :);
full_csd  = Crow * X_ar;

% ---- SPLIT PIPELINE: EXACTLY what sim_inject_power.m now computes ----
% This must mirror the production expression term for term. The previous
% version of this script average-referenced the injection as a full matrix and
% so tested correct algebra that the production code did not implement: that
% code added the RAW injection rows to average-referenced real data, which is
% right only if the leadfield gains are zero-mean across channels -- asserted
% in a comment, never measured. The gains are now average-referenced explicitly
% (gFa, gTa) and projected once, and this check uses the same vectors.
Xr             = X_real - mean(X_real, 1);       % avg-ref of the real data only
split_real_avg = Xr(iFwd, :);
split_real_csd = Crow * Xr;

gFa = gF - mean(gF);                              % as sim_inject_power.m
gTa = gT - mean(gT);
pF_avg = sF_sc * gFa(iFwd);
pT_avg = sT_sc * gTa(iFwd);
pF_csd = sF_sc * (Crow * gFa);
pT_csd = sT_sc * (Crow * gTa);

split_avg = split_real_avg + pF_avg * sA + pT_avg * s2;
split_csd = split_real_csd + pF_csd * sA + pT_csd * s2;

% ---- how far the gains are from zero-mean: the assumption that was never tested ----
gmean_rel = max(abs([mean(gF) mean(gT)])) / max(max(abs([gF(:); gT(:)])), realmin);
fprintf('\ngain vectors: |mean| / max|g| = %.3e\n', gmean_rel);
fprintf('   (if this were ~0 the previous version was also correct; it is reported, not assumed)\n');

% ---- COMPARE ----
tol = 1e-9;
rng_avg = max(abs(full_avg(:)));
rng_csd = max(abs(full_csd(:)));

err_avg = max(abs(full_avg(:) - split_avg(:))) / rng_avg;
err_csd = max(abs(full_csd(:) - split_csd(:))) / rng_csd;

fprintf('\n--- Linearity check results ---\n');
fprintf('  (comparing against the expression sim_inject_power.m actually evaluates)\n');
fprintf('  Average-ref ROI:  max relative error = %.3e  (tol = %.0e)   %s\n', ...
    err_avg, tol, local_verdict(err_avg, tol));
fprintf('  CSD ROI:          max relative error = %.3e  (tol = %.0e)   %s\n', ...
    err_csd, tol, local_verdict(err_csd, tol));

if err_avg < tol && err_csd < tol
    fprintf('\n*** PASSED: linearity shortcut is algebraically exact. ***\n\n');
    f_out = fullfile(cfg.out_dir, 'sim_inject_linearity_check.md');
    fid = fopen(f_out, 'w');
    fprintf(fid, '# sim_inject_linearity_check\n\n');
    fprintf(fid, 'Subject: %s | Run: %s\n\n', subj, key);
    fprintf(fid, '| Check | Max relative error | Tolerance | Verdict |\n');
    fprintf(fid, '|-------|-------------------|-----------|--------|\n');
    fprintf(fid, '| Average-ref ROI | %.3e | %.0e | PASSED |\n', err_avg, tol);
    fprintf(fid, '| CSD ROI         | %.3e | %.0e | PASSED |\n', err_csd, tol);
    fprintf(fid, '\nGain vectors |mean|/max|g| = %.3e. The split path average-references the\n', gmean_rel);
    fprintf(fid, 'two gain vectors explicitly, so this number is reported rather than relied on.\n');
    fprintf(fid, 'The split path mirrors sim_inject_power.m term for term (rank-2 outer products,\n');
    fprintf(fid, 'no full nchan x nsamp injection matrix).\n');
    fprintf(fid, '\nResult: PASSED\n');
    fclose(fid);
    fprintf('Results saved to %s\n', f_out);
else
    error('FAILED: linearity check did not pass. DO NOT run sim_inject_power with the shortcut.');
end

function idx = local_idx(have, want, where)
    have = reshape(cellstr(have), 1, []);
    idx = zeros(1, numel(want));
    for i = 1:numel(want)
        j = find(strcmpi(have, want{i}), 1);
        if isempty(j), error('channel "%s" not present in the %s.', want{i}, where); end
        idx(i) = j;
    end
end

function y = local_delay(x, tau)
    if tau == 0, y = x; return; end
    n = numel(x);
    X = fft(x);
    k = 0:n-1;
    f = k; f(k > n/2) = k(k > n/2) - n;
    y = real(ifft(X .* exp(-2i*pi*f*tau/n)));
    s = std(y); if s > 0, y = y / s; end
end

function v = local_verdict(err, tol)
    if err < tol, v = 'PASSED'; else, v = 'FAILED'; end
end
