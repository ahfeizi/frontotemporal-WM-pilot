% sim_inject_power.m -- Dataset 2 injection power curve.
%
% Injects a forward-projected lagged frontal(F3)-temporal(T7) source pair into
% the REAL preprocessed high-load trials at five coupling strengths and asks how
% often the battery notices. N = 22, each subject's real trial counts, each
% subject's own runs pooled as s03_connectivity.m pools them.
%
% The two arms get the SAME real trials, the SAME two sources and the SAME
% amplitude; only the coupling differs, and source 2's mix is variance
% preserving, so the injected power is identical in both arms at every coupling
% level. kappa = 0 is therefore the false-positive row: the arms then differ
% only by two independent noise realisations of source 2. Because both arms
% share one real recording, the paired comparison is as sensitive as the design
% allows, which is what a power analysis wants.
%
% This is a power analysis WITHIN high load, not a re-run of the high-low
% contrast. The beta band receives no injection and is the specificity check.
%
% Injection amplitude is cfg.inject_frac of each subject's OWN observed ROI
% amplitude in the target band, measured from that subject's real data.
%
% RESUME: each subject is appended to sim_inject_contrasts.csv and then listed
% in sim_inject_done.txt. Re-running skips completed subjects and strips partial
% rows. Inference runs in sim_inject_curve.m, once all subjects are present.
%
% Run:  cd <P1>/revision/sim; sim_inject_power
% Then: sim_inject_curve
%
% LINEARITY SHORTCUT (approved 2026-10-30).
% Because average-reference, CSD, filtfilt and Hilbert/wavelet are all linear,
% the analytic signal of (real + injection) equals the analytic signal of the
% real data plus that of the injection alone.  The real part is therefore pre-
% computed ONCE per subject per run, before the repetition loop, and cached as:
%
%   rdata(k).avg_roi   -- Ch_roi x T, average-ref ROI rows of the real data
%   rdata(k).csd_roi   -- Ch_roi x T, CSD ROI rows of the real data
% and, because the injection is rank two, its two sensor patterns are
% average-referenced and CSD-projected once per subject as well, so no full
% nchan x nsamp injection matrix is ever formed.
%
% Every repetition then only processes the injected source signals (tiny vectors)
% and adds the precomputed real part, cutting this stage by roughly x6.
%
% Algebraic identity verified by sim_inject_linearity_check.m (run once, must PASS
% before sim_inject_power).  Shape: pick one subject, one run, one kappa; compare
% the full pipeline against the split pipeline for a single rep.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = sim_inject_config();
repo = cfg.repo_root;
addpath(repo);                              % config.m, util_bandpass.m, util_chanidx.m, util_conn_pairs.m
addpath(cfg.shared_sim_dir);                % sim_util_path.m, sim_util_csd.m, sim_util_stats.m
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end
if ~isfolder(cfg.cache_dir), mkdir(cfg.cache_dir); end

sim_util_path(cfg, {'ft_defaults','ft_prepare_leadfield','ft_prepare_sourcemodel', ...
                    'ft_scalpcurrentdensity','ft_convert_units','ft_datatype_sens'});

f_main = fullfile(cfg.out_dir, 'sim_inject_contrasts.csv');
f_done = fullfile(cfg.out_dir, 'sim_inject_done.txt');
hdr = ['subject,kappa,rep,band,measure,decomposition,reference,contrast,' ...
       'roimean_armA,roimean_armB,ntrials_used,obs_roi_sd,inject_scale_frontal,' ...
       'inject_scale_temporal,tau_samples,seed'];
done = local_resume(f_done, f_main, hdr);

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
fprintf('%d subjects in the manifest\n', numel(subjects));

win  = round(cfg.maint_window * cfg.srate);
wlen = win(2) - win(1) + 1;
tb   = cfg.target_band;
band = cfg.bands.(tb);
fc   = mean(band);
tau  = (cfg.phase_offset_deg/360) * cfg.srate / fc;      % samples
nK   = numel(cfg.kappas);
nRep = cfg.n_rep_hilbert;
nWav = min(cfg.n_rep_wavelet, nRep);

% Built once, from the first subject actually processed, and reused for the rest
% (the montage is asserted identical for every subsequent subject).
F = [];        % forward model
Crow = [];     % ROI rows of the cached CSD operator
iFwd = [];     % ROI rows of the full montage
gF = []; gT = [];   % sensor patterns of the two injection sites

for su = 1:numel(subjects)
    subj = char(subjects(su));
    if ismember(subj, done)
        fprintf('%s already complete, skipping\n', subj);
        continue
    end
    t0 = tic;
    keys = manifest.key(manifest.subject == subjects(su));

    % ---------------- load this subject's runs ----------------
    runs = struct('srate',{},'locks',{},'full',{});
    labF = {}; labT = {}; roiLabs = {};
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if k == 1
            have = {EEG.chanlocs.labels};
            labF = cfg.roi_frontal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_frontal));
            labT = cfg.roi_temporal(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_temporal));
            if numel(labF) < numel(cfg.roi_frontal) || numel(labT) < numel(cfg.roi_temporal)
                error('%s: only %d/%d frontal and %d/%d temporal ROI channels present.', ...
                    subj, numel(labF), numel(cfg.roi_frontal), numel(labT), numel(cfg.roi_temporal));
            end
            roiLabs = [labF labT];
            % forward model and CSD operator, once for the whole dataset
            ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
            if isempty(F)
                F = sim_inject_forward(cfg, ft_raw.elec, ft_raw.label, 'ds004117');
                [Cfull, csd_labels] = sim_util_csd(cfg, F.elec, F.labels, cfg.srate, 'ds004117');
                iCsd = local_idx(csd_labels, roiLabs, 'CSD output');
                Crow = Cfull(iCsd, :);
                iFwd = local_idx(F.labels, roiLabs, 'forward montage');
                gF = F.sites.frontal.gain;
                gT = F.sites.temporal.gain;
            end
            if ~isequal(reshape(cellstr(ft_raw.label),1,[]), F.labels)
                error('%s: montage differs from the one the forward model was built for.', subj);
            end
        end
        sel = trials([trials.good] & strcmp({trials.load}, cfg.load_condition));
        runs(k).srate = EEG.srate; %#ok<SAGROW>
        runs(k).locks = [sel.lock_sample]; %#ok<SAGROW>
        runs(k).full  = double(EEG.data); %#ok<SAGROW>
        clear EEG
    end
    if isempty(iFwd) || isempty(Crow) || isempty(gF) || isempty(gT)
        error('%s: the forward model / CSD operator setup did not run.', subj);
    end
    ntr_total = sum(cellfun(@numel, {runs.locks}));
    if ntr_total < cfg.min_trials
        warning('%s: only %d good %s-load trials, skipping subject.', subj, ntr_total, cfg.load_condition);
        continue
    end

    % ---------------- observed ROI amplitude in the target band ----------------
    % Measured on this subject's own real data, average reference, ROI channels.
    sds = [];
    for k = 1:numel(runs)
        Xr = runs(k).full - mean(runs(k).full, 1);
        Zb = util_bandpass(Xr(iFwd,:), runs(k).srate, band, cfg.filt_order_factor);
        sds = [sds std(Zb, 0, 2)']; %#ok<AGROW>
    end
    obs_roi_sd = mean(sds);
    sF = cfg.inject_frac * obs_roi_sd / mean(abs(gF(iFwd)));
    sT = cfg.inject_frac * obs_roi_sd / mean(abs(gT(iFwd)));
    fprintf('%s: %d %s-load trials over %d run(s), observed ROI %s SD %.4g\n', ...
        subj, ntr_total, cfg.load_condition, numel(runs), tb, obs_roi_sd);

    nF = numel(labF); nT = numel(labT);
    fr = 1:nF; te = nF + (1:nT);

    % ---- LINEARITY SHORTCUT: precompute real-data ROI rows (once per subject) ----
    % rdata(k).avg_roi  -- average-ref ROI rows of the real data (iFwd rows of X_avgref)
    % rdata(k).csd_roi  -- CSD ROI rows of the real data (Crow * X_avgref)
    % These are added to the injection's own ROI rows inside local_arm.
    rdata = struct('avg_roi',{},'csd_roi',{});
    for k = 1:numel(runs)
        Xr = runs(k).full - mean(runs(k).full, 1);   % average reference
        rdata(k).avg_roi = Xr(iFwd, :);              %#ok<SAGROW>
        rdata(k).csd_roi = Crow * Xr;                %#ok<SAGROW>
    end

    % ---- the injection is RANK TWO, so project its two gain vectors once ----
    % inj(t) = sF*gF*sA(t) + sT*gT*s2(t), i.e. two fixed sensor patterns times
    % two scalar time courses. Everything linear applied to it can therefore be
    % applied to gF and gT alone, ONCE per subject, instead of to a full
    % nchan x nsamp matrix once per arm. At 69 channels and 2,103,619 samples
    % that matrix is 1.16 GB, built and discarded six times per repetition.
    %
    % The average reference is taken HERE, on the gains, rather than assumed
    % away. The earlier version of this file added the raw injection rows to
    % already-average-referenced real data, which is only correct if the
    % leadfield gains happen to be zero-mean across channels -- an assumption
    % that was stated in a comment and never tested. Subtracting the mean from
    % the gain vectors makes it true by construction and costs two scalars.
    gFa = gF - mean(gF);                 % average-referenced sensor patterns
    gTa = gT - mean(gT);
    pF_avg = sF * gFa(iFwd);             % nroi x 1, average-reference ROI rows
    pT_avg = sT * gTa(iFwd);
    pF_csd = sF * (Crow * gFa);          % nroi x 1, CSD ROI rows
    pT_csd = sT * (Crow * gTa);

    % ---------------- repetitions ----------------
    buf = cell(1, nRep);
    for rep = 1:nRep
        buf{rep} = local_one_rep(cfg, runs, rdata, F, pF_avg, pT_avg, pF_csd, pT_csd, sF, sT, ...
            band, tau, win, wlen, fr, te, subj, su, rep, rep <= nWav);
        if mod(rep, 25) == 0
            fprintf('   %s rep %d/%d (%.1f min)\n', subj, rep, nRep, toc(t0)/60);
        end
    end

    rows = [buf{:}];
    local_append(f_main, rows);
    local_append(f_done, {sprintf('%s\n', subj)});
    fprintf('%s complete: %d rows (%.1f min)\n', subj, numel(rows), toc(t0)/60);
    clear runs rdata
end

fprintf('\nAll subjects complete. %s\nNow run sim_inject_curve.m\n', f_main);

% ======================= local functions =======================
function rows = local_one_rep(cfg, runs, rdata, F, pF_avg, pT_avg, pF_csd, pT_csd, sF, sT, ...
        band, tau, win, wlen, fr, te, subj, su, rep, do_wavelet)
    % One repetition: arm A once, then arm B at each coupling level, all sharing
    % the same source-1 realisation so the arms differ only as intended.
    %
    % LINEARITY SHORTCUT: only injection source signals (1 x T vectors) are
    % bandpass-filtered here.  The real-data ROI rows (rdata) were precomputed
    % before this loop and are simply added to the injected contribution.
    seed = cfg.seed_base + 1e6*su + 1e3*rep;
    rng(seed, 'twister');

    nrun = numel(runs);
    sA = cell(1,nrun); nA2 = cell(1,nrun); nB2 = cell(1,nrun); sAd = cell(1,nrun);
    for k = 1:nrun
        n = size(runs(k).full, 2);
        sA{k}  = local_bandnoise(n, runs(k).srate, band, cfg.filt_order_factor);
        nA2{k} = local_bandnoise(n, runs(k).srate, band, cfg.filt_order_factor);
        nB2{k} = local_bandnoise(n, runs(k).srate, band, cfg.filt_order_factor);
        sAd{k} = local_delay(sA{k}, tau);
    end

    % arm A: coupling 0, its own source-2 noise
    armA = local_arm(cfg, runs, rdata, pF_avg, pT_avg, pF_csd, pT_csd, sA, nA2, [], 0, win, wlen, fr, te, do_wavelet);

    rows = {};
    for ik = 1:numel(cfg.kappas)
        kap = cfg.kappas(ik);
        armB = local_arm(cfg, runs, rdata, pF_avg, pT_avg, pF_csd, pT_csd, sA, nB2, sAd, kap, ...
            win, wlen, fr, te, do_wavelet);
        if do_wavelet, decomps = {'hilbert','wavelet'}; else, decomps = {'hilbert'}; end
        for r = 1:numel(cfg.refs)
            rf = cfg.refs{r};
            for dd = 1:numel(decomps)
                dc = decomps{dd};
                for b = 1:numel(cfg.bands_used)
                    bn = cfg.bands_used{b};
                    for m = 1:numel(cfg.measures)
                        me = cfg.measures{m};
                        a  = mean(armA.(rf).(dc).(bn).(me), 'omitnan');
                        bb = mean(armB.(rf).(dc).(bn).(me), 'omitnan');
                        rows{end+1} = sprintf(['%s,%.17g,%d,%s,%s,%s,%s,%.17g,%.17g,%.17g,' ...
                            '%d,%.17g,%.17g,%.17g,%.17g,%d\n'], ...
                            subj, kap, rep, bn, me, dc, rf, bb - a, a, bb, ...
                            armB.ntrials_used, cfg.inject_frac, sF, sT, tau, seed); %#ok<AGROW>
                    end
                end
            end
        end
    end
end

function V = local_arm(cfg, runs, rdata, pF_avg, pT_avg, pF_csd, pT_csd, sA, n2, sAd, kap, ...
        win, wlen, fr, te, do_wavelet)
    % Build this arm's ROI-row data for every run, then run the battery.
    %
    % LINEARITY SHORTCUT, in the form the rank-2 structure allows. Average
    % reference, the CSD operator, filtfilt and hilbert are all linear, so
    %   avg_ref(real + inj) = avg_ref(real) + avg_ref(inj)
    %   CSD    (real + inj) = CSD(real)     + CSD(inj)
    % with avg_ref(real) and CSD(real) precomputed once per subject (rdata), and
    % the injection's two sensor patterns already average-referenced and
    % CSD-projected once per subject (pF_*, pT_*).
    %
    % What is left here is two outer products per reference: an nroi x 1 vector
    % times a 1 x nsamp time course. The full nchan x nsamp injection matrix is
    % never formed, so this neither allocates 1.16 GB per arm nor multiplies the
    % CSD operator through it -- about 11 x 2.1M flops per term in place of
    % 69 x 2.1M element writes plus a 1.6 GFLOP operator multiply.
    arm = struct('srate',{},'locks',{},'average',{},'csd',{});
    for k = 1:numel(runs)
        if isempty(sAd)
            s2 = n2{k};                                        % coupling 0
        else
            s2 = kap * sAd{k} + sqrt(max(1-kap^2,0)) * n2{k};  % variance preserving
        end
        arm(k).srate   = runs(k).srate; %#ok<SAGROW>
        arm(k).locks   = runs(k).locks; %#ok<SAGROW>
        arm(k).average = rdata(k).avg_roi + pF_avg * sA{k} + pT_avg * s2; %#ok<SAGROW>
        arm(k).csd     = rdata(k).csd_roi + pF_csd * sA{k} + pT_csd * s2; %#ok<SAGROW>
    end
    V = sim_inject_battery(cfg, arm, win, wlen, fr, te, do_wavelet);
end

function s = local_bandnoise(n, srate, band, order_factor)
    s = util_bandpass(randn(1,n), srate, band, order_factor);
    sd = std(s);
    if sd > 0, s = s / sd; end
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

function idx = local_idx(have, want, where)
    have = reshape(cellstr(have), 1, []);
    idx = zeros(1, numel(want));
    for i = 1:numel(want)
        j = find(strcmpi(have, want{i}), 1);
        if isempty(j), error('channel "%s" not present in the %s.', want{i}, where); end
        idx(i) = j;
    end
end

function local_append(f, buf)
    fid = fopen(f, 'a');
    if fid < 0, error('cannot open %s for appending', f); end
    fprintf(fid, '%s', buf{:});
    fclose(fid);
end

function done = local_resume(f_done, f_main, hdr)
    done = {};
    if isfile(f_done)
        done = strtrim(splitlines(fileread(f_done)));
        done = reshape(done(~cellfun(@isempty, done)), 1, []);
    end
    if ~isfile(f_main)
        if ~isempty(done)
            error('%s is missing but %s lists completed subjects. Delete %s to restart.', ...
                f_main, f_done, f_done);
        end
        fid = fopen(f_main,'w'); fprintf(fid,'%s\n', hdr); fclose(fid);
        return
    end
    L = splitlines(fileread(f_main));
    L = L(~cellfun(@isempty, L));
    body = L(2:end);
    subj = regexp(body, '^[^,]*', 'match', 'once');
    keep = ismember(subj, done);
    if all(keep)
        if ~isempty(done), fprintf('Resume: %d subject(s) already complete.\n', numel(done)); end
        return
    end
    fprintf('Resume: removing %d row(s) of an unfinished subject from %s\n', sum(~keep), f_main);
    fid = fopen(f_main,'w'); fprintf(fid,'%s\n', L{1}, body{keep}); fclose(fid);
end
