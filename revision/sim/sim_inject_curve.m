% sim_inject_curve.m -- inference and aggregation for the injection power curve.
%
% Reads sim_inject_contrasts.csv (one row per subject x kappa x rep x cell),
% forms the per-subject contrast vector of every cell, runs the EXISTING
% inference on each (util_contrast_stats.contrast_stats through
% sim_stats_runner.py, so the sign-flip test is the exact enumeration of all
% 2^22 flips at N = 22), and aggregates into the power curve:
%
%   detection probability = P(sign-flip p < cfg.alpha) per coupling level x
%   measure x decomposition x reference, with the kappa = 0 row as the
%   false-positive rate and P(TOST at +/-0.4 passes) reported there too.
%
% The beta band carries no injection, so its rows are the specificity check.
%
% Outputs (revision/sim/results/):
%   sim_inject_percell.csv   one row per kappa x rep x cell, full inference
%   sim_inject_curve.csv     the aggregated curve
%   sim_inject_power.md      at most 20 lines
%
% Run:  cd <P1>/revision/sim; sim_inject_curve

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = sim_inject_config();
addpath(cfg.shared_sim_dir);                % sim_util_stats.m

f_in = fullfile(cfg.out_dir, 'sim_inject_contrasts.csv');
if ~isfile(f_in), error('%s not found. Run sim_inject_power.m first.', f_in); end
T = readtable(f_in, 'TextType', 'string');
fprintf('read %d rows from %s\n', height(T), f_in);

subjects = unique(T.subject, 'stable');
nsub = numel(subjects);
fprintf('%d subjects\n', nsub);

% ---------------- one cell per kappa x rep x band x measure x decomp x ref ----
G = [string(T.kappa), string(T.rep), T.band, T.measure, T.decomposition, T.reference];
[u, ~, gi] = unique(G, 'rows');
nG = size(u,1);
fprintf('%d cells to infer at N = %d\n', nG, nsub);

keys = cell(1, nG);
X = nan(nG, nsub);
ntr = zeros(1, nG);
for g = 1:nG
    k = find(gi == g);
    [tf, loc] = ismember(T.subject(k), subjects);
    if ~all(tf), error('cell %d references an unknown subject.', g); end
    X(g, loc(tf)) = T.contrast(k(tf));
    ntr(g) = max(T.ntrials_used(k));
    keys{g} = sprintf('%s|%s|%s|%s|%s|%s', u(g,1), u(g,2), u(g,3), u(g,4), u(g,5), u(g,6));
end

R = sim_util_stats(cfg, keys, X, 'inject');

% ---------------- per-cell output ----------------
f_cell = fullfile(cfg.out_dir, 'sim_inject_percell.csv');
fid = fopen(f_cell,'w');
fprintf(fid, ['kappa,rep,band,measure,decomposition,reference,n,dz,dz_ci_lo,dz_ci_hi,' ...
    't,df,p,p_perm,perm_method,bf10,bf01,p_tost,mean_contrast,sd_contrast,' ...
    'ntrials_used_max,inject_frac,seed_base\n']);
for i = 1:nG
    p = strsplit(keys{i}, '|');
    fprintf(fid, ['%s,%s,%s,%s,%s,%s,%d,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%s,' ...
        '%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%d\n'], ...
        p{1}, p{2}, p{3}, p{4}, p{5}, p{6}, R.n(i), R.dz(i), R.dz_ci_lo(i), R.dz_ci_hi(i), ...
        R.t(i), R.df(i), R.p(i), R.p_perm(i), char(R.perm_method(i)), ...
        R.bf10(i), R.bf01(i), R.p_tost(i), R.mean(i), R.sd(i), ntr(i), ...
        cfg.inject_frac, cfg.seed_base);
end
fclose(fid);
fprintf('wrote %s (%d rows)\n', f_cell, nG);

% ---------------- aggregate into the curve ----------------
kap = str2double(u(:,1));   % string -> number, unambiguous for a string array
A = [string(kap), u(:,3), u(:,4), u(:,5), u(:,6)];
[ua, ~, ai] = unique(A, 'rows');
nA = size(ua,1);
curve = struct('kappa',{},'band',{},'measure',{},'decomposition',{},'reference',{}, ...
               'n_rep',{},'p_detect',{},'p_bf',{},'p_tost',{},'dz_mean',{});
rows = cell(1, nA);
for a = 1:nA
    k = (ai == a);
    rec = struct();
    rec.kappa = str2double(ua(a,1));
    rec.band = ua(a,2); rec.measure = ua(a,3);
    rec.decomposition = ua(a,4); rec.reference = ua(a,5);
    rec.n_rep = sum(k);
    rec.p_detect = mean(R.p_perm(k) < cfg.alpha);
    rec.p_detect_par = mean(R.p(k) < cfg.alpha);
    rec.p_bf = mean(R.bf10(k) > cfg.bf_thresh);
    rec.p_tost = mean(R.p_tost(k) < cfg.alpha);
    rec.dz_mean = mean(R.dz(k), 'omitnan');
    rec.dz_sd = std(R.dz(k), 'omitnan');
    curve(end+1) = struct('kappa',rec.kappa,'band',rec.band,'measure',rec.measure, ...
        'decomposition',rec.decomposition,'reference',rec.reference,'n_rep',rec.n_rep, ...
        'p_detect',rec.p_detect,'p_bf',rec.p_bf,'p_tost',rec.p_tost,'dz_mean',rec.dz_mean); %#ok<SAGROW>
    rows{a} = sprintf('%.17g,%s,%s,%s,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.17g,%.17g\n', ...
        rec.kappa, rec.band, rec.measure, rec.decomposition, rec.reference, ...
        rec.n_rep, nsub, rec.p_detect, rec.p_detect_par, rec.p_bf, rec.p_tost, ...
        rec.dz_mean, rec.dz_sd);
end

f_curve = fullfile(cfg.out_dir, 'sim_inject_curve.csv');
fid = fopen(f_curve,'w');
fprintf(fid, ['kappa,band,measure,decomposition,reference,n_rep,n_subjects,' ...
    'p_detect,p_detect_parametric,p_bf10_gt3,p_tost_pass,dz_mean,dz_sd\n']);
fprintf(fid, '%s', rows{:});
fclose(fid);
fprintf('wrote %s (%d rows)\n', f_curve, nA);

% ---------------- 20-line summary ----------------
tb = cfg.target_band;
L = {};
L{end+1} = '# Dataset 2 injection power curve (sim_inject_curve.m)';
L{end+1} = sprintf(['Lagged %s-%s pair injected into the real %s-load trials, N = %d, real trial ' ...
    'counts, %d reps per level.'], cfg.site_frontal, cfg.site_temporal, cfg.load_condition, ...
    nsub, cfg.n_rep_hilbert);
L{end+1} = sprintf('Detection = P(sign-flip p < %.2f). Band %s carries the injection; beta carries none.', ...
    cfg.alpha, tb);
for r = 1:numel(cfg.refs)
    L{end+1} = sprintf('## %s reference, Hilbert, band %s', cfg.refs{r}, tb);
    for m = 1:numel(cfg.measures)
        s = '';
        for ik = 1:numel(cfg.kappas)
            v = local_get(curve, cfg.kappas(ik), tb, cfg.measures{m}, 'hilbert', cfg.refs{r}, 'p_detect');
            s = [s sprintf('k=%.3g %s  ', cfg.kappas(ik), local_pc(v))]; %#ok<AGROW>
        end
        L{end+1} = sprintf('- %-5s %s', cfg.measures{m}, strtrim(s));
    end
end
L{end+1} = sprintf('## kappa = 0 (no coupling difference): false positives and equivalence');
for r = 1:numel(cfg.refs)
    s = '';
    for m = 1:numel(cfg.measures)
        fp = local_get(curve, 0, tb, cfg.measures{m}, 'hilbert', cfg.refs{r}, 'p_detect');
        tp = local_get(curve, 0, tb, cfg.measures{m}, 'hilbert', cfg.refs{r}, 'p_tost');
        s = [s sprintf('%s FP %s / TOST %s  ', cfg.measures{m}, local_pc(fp), local_pc(tp))]; %#ok<AGROW>
    end
    L{end+1} = sprintf('- %s: %s', cfg.refs{r}, strtrim(s));
end
L{end+1} = sprintf('Specificity: beta PLV detection at the strongest level, average ref = %s.', ...
    local_pc(local_get(curve, cfg.kappas(end), 'beta', 'PLV', 'hilbert', 'average', 'p_detect')));
L{end+1} = sprintf('Observed, for comparison: theta PLV high-low dz %+.3f, BF01 %.2f.', ...
    cfg.obs.theta_plv_dz, cfg.obs.theta_plv_bf01);
L{end+1} = 'This is a power analysis within high load, not a re-run of the high-low contrast.';
L{end+1} = sprintf('Per-cell inference %s; full curve %s.', f_cell, f_curve);

f_md = fullfile(cfg.out_dir, 'sim_inject_power.md');
fid = fopen(f_md,'w'); fprintf(fid,'%s\n', L{:}); fclose(fid);
fprintf('wrote %s (%d lines)\n', f_md, numel(L));

% ======================= local functions =======================
function v = local_get(curve, kap, bn, me, dc, rf, field)
    k = find(abs([curve.kappa] - kap) < 1e-12 & ...
             ismember([curve.band]', string(bn)) & ...
             ismember([curve.measure]', string(me)) & ...
             ismember([curve.decomposition]', string(dc)) & ...
             ismember([curve.reference]', string(rf)), 1);
    if isempty(k), v = NaN; else, v = curve(k).(field); end
end

function s = local_pc(v)
    if isnan(v), s = '--'; else, s = sprintf('%.0f%%', 100*v); end
end
