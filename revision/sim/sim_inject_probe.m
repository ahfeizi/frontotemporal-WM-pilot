% sim_inject_probe.m -- measure the injection power curve's cost, then project.
%
% Runs cfg.probe_reps repetitions for the first cfg.probe_subjects subjects and
% times each stage: loading a subject, building the forward model and CSD
% operator (once), one arm's battery on the Hilbert path, one arm's battery with
% the wavelet path, and the batched inference. Then projects the whole stage
% from the counts in sim_inject_config.m.
%
% Writes revision/sim/results/sim_inject_timing.md and prints the same.
%
% Run:  cd <P1>/revision/sim; sim_inject_probe

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = sim_inject_config();
addpath(cfg.repo_root);
addpath(cfg.shared_sim_dir);
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end
if ~isfolder(cfg.cache_dir), mkdir(cfg.cache_dir); end

sim_util_path(cfg, {'ft_defaults','ft_prepare_leadfield','ft_prepare_sourcemodel', ...
                    'ft_scalpcurrentdensity','ft_convert_units','ft_datatype_sens'});

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
nsub_all = numel(subjects);
nprobe = min(cfg.probe_subjects, nsub_all);

win  = round(cfg.maint_window * cfg.srate);
wlen = win(2) - win(1) + 1;
band = cfg.bands.(cfg.target_band);
tau  = (cfg.phase_offset_deg/360) * cfg.srate / mean(band);

t_load = []; t_hil = []; t_wav = []; t_prep = NaN;
nsamp_total = 0; ntr_total = 0;
F = []; Crow = []; iFwd = []; gF = []; gT = [];

for su = 1:nprobe
    subj = char(subjects(su));
    keys = manifest.key(manifest.subject == subjects(su));
    tl = tic;
    runs = struct('srate',{},'locks',{},'full',{});
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if k == 1 && isempty(F)
            have = {EEG.chanlocs.labels};
            labF = cfg.roi_frontal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_frontal));
            labT = cfg.roi_temporal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_temporal));
            roiLabs = [labF labT];
            ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
            tp = tic;
            F = sim_inject_forward(cfg, ft_raw.elec, ft_raw.label, 'ds004117');
            [Cfull, csd_labels] = sim_util_csd(cfg, F.elec, F.labels, cfg.srate, 'ds004117');
            t_prep = toc(tp);
            iCsd = local_idx(csd_labels, roiLabs);
            Crow = Cfull(iCsd, :);
            iFwd = local_idx(F.labels, roiLabs);
            gF = F.sites.frontal.gain; gT = F.sites.temporal.gain;
            nF = numel(labF); nT = numel(labT);
            fr = 1:nF; te = nF + (1:nT);
        end
        sel = trials([trials.good] & strcmp({trials.load}, cfg.load_condition));
        runs(k).srate = EEG.srate; %#ok<SAGROW>
        runs(k).locks = [sel.lock_sample]; %#ok<SAGROW>
        runs(k).full  = double(EEG.data); %#ok<SAGROW>
        nsamp_total = nsamp_total + size(EEG.data,2);
        clear EEG
    end
    t_load(end+1) = toc(tl); %#ok<SAGROW>
    ntr_total = ntr_total + sum(cellfun(@numel, {runs.locks}));

    for rep = 1:cfg.probe_reps
        rng(cfg.seed_base + su*1000 + rep, 'twister');
        nrun = numel(runs);
        sA = cell(1,nrun); n2 = cell(1,nrun);
        for k = 1:nrun
            n = size(runs(k).full,2);
            sA{k} = local_bn(n, runs(k).srate, band, cfg.filt_order_factor);
            n2{k} = local_bn(n, runs(k).srate, band, cfg.filt_order_factor);
        end
        arm = struct('srate',{},'locks',{},'average',{},'csd',{});
        for k = 1:nrun
            X = runs(k).full + gF * sA{k} + gT * n2{k};
            X = X - mean(X,1);
            arm(k).srate = runs(k).srate; %#ok<SAGROW>
            arm(k).locks = runs(k).locks; %#ok<SAGROW>
            arm(k).average = X(iFwd,:); %#ok<SAGROW>
            arm(k).csd = Crow * X; %#ok<SAGROW>
        end
        th = tic; sim_inject_battery(cfg, arm, win, wlen, fr, te, false); t_hil(end+1) = toc(th); %#ok<SAGROW>
        tw = tic; sim_inject_battery(cfg, arm, win, wlen, fr, te, true);  t_wav(end+1) = toc(tw); %#ok<SAGROW>
    end
    fprintf('  %s timed (%d rep(s))\n', subj, cfg.probe_reps);
    clear runs
end

% ---- inference, measured on a real batched call ----
nb = 12;
kk = arrayfun(@(k) sprintf('probe_%d', k), 1:nb, 'UniformOutput', false);
rng(1,'twister'); Xb = 0.3 + randn(nb, nsub_all);
ti = tic; sim_util_stats(cfg, kk, Xb, 'inject_probe'); t_inf = toc(ti);
inf_per_cell = t_inf / nb;

% ---------------- projection ----------------
nK = numel(cfg.kappas);
nRep = cfg.n_rep_hilbert; nWav = min(cfg.n_rep_wavelet, nRep);
arms_per_rep = 1 + nK;                      % arm A once, arm B per coupling level
sim_sec = nsub_all * ( (nRep - nWav) * arms_per_rep * mean(t_hil) + ...
                       nWav * arms_per_rep * mean(t_wav) );
load_sec = nsub_all * mean(t_load);
ncell = nK * (nRep + nWav) * numel(cfg.bands_used) * numel(cfg.measures) * numel(cfg.refs);
inf_sec = ncell * inf_per_cell;

L = {};
L{end+1} = '# Injection power curve: run-time probe (sim_inject_probe.m)';
L{end+1} = sprintf('Measured on this machine, %s, over %d subject(s) x %d rep(s).', ...
    datestr(now,'yyyy-mm-dd'), nprobe, cfg.probe_reps);
L{end+1} = sprintf('- forward model + CSD operator: %.1f s, once (cached afterwards)', t_prep);
L{end+1} = sprintf('- loading one subject: %.1f s (%d samples, %d %s-load trials over the probe set)', ...
    mean(t_load), nsamp_total, ntr_total, cfg.load_condition);
L{end+1} = sprintf('- one arm, Hilbert path: %.2f s; one arm, Hilbert+wavelet: %.2f s', ...
    mean(t_hil), mean(t_wav));
L{end+1} = sprintf('- inference: %.3f s per cell at N = %d', inf_per_cell, nsub_all);
L{end+1} = sprintf('Counts: %d coupling levels, %d reps (%d also on the wavelet path), %d arms per rep.', ...
    nK, nRep, nWav, arms_per_rep);
L{end+1} = sprintf('**Projection: %.1f h simulation + %.1f h inference + %.1f min loading = %.1f h total.**', ...
    sim_sec/3600, inf_sec/3600, load_sec/60, (sim_sec+inf_sec+load_sec)/3600);
L{end+1} = sprintf('%d inference cells in all; sim_stats_runner.py parallelises them internally.', ncell);
L{end+1} = 'The simulation part is a plain subject loop and is not parallelised, because each';
L{end+1} = 'subject holds its whole recording in memory. Run it overnight, or lower';
L{end+1} = 'cfg.n_rep_hilbert / cfg.n_rep_wavelet, which scale the total linearly.';
L{end+1} = 'The linearity shortcut noted at the top of sim_inject_power.m would cut the';
L{end+1} = 'simulation part by roughly an order of magnitude if you want it.';

f_md = fullfile(cfg.out_dir, 'sim_inject_timing.md');
fid = fopen(f_md,'w'); fprintf(fid,'%s\n', L{:}); fclose(fid);
fprintf('\n%s\n', L{:});
fprintf('\nwrote %s\n', f_md);

% ======================= local functions =======================
function idx = local_idx(have, want)
    have = reshape(cellstr(have), 1, []);
    idx = zeros(1, numel(want));
    for i = 1:numel(want)
        j = find(strcmpi(have, want{i}), 1);
        if isempty(j), error('channel "%s" not present.', want{i}); end
        idx(i) = j;
    end
end

function s = local_bn(n, srate, band, order_factor)
    s = util_bandpass(randn(1,n), srate, band, order_factor);
    sd = std(s); if sd > 0, s = s / sd; end
end
