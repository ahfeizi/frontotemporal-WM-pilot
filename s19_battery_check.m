% s19_battery_check.m -- Regression check of s18_battery.m output against
% the existing pipeline results. Reads battery_persubject.csv and
% battery_phasediff.csv from cfg.deriv_dir and prints PASS/FAIL per check.
%
%   T1  Hilbert / average-ref ROI-average PLV and dwPLI per subject x load
%       == conn_results.csv           (theta, beta;          s03)
%   T2  same                          == conn_results_betasplit.csv
%                                        (beta_low, beta_high; s09)
%   T3  Wavelet / average-ref ROI-average PLV and dwPLI per subject x load
%       == conn_results_wavelet.csv   (bands present in that file; s16 did
%                                        not compute theta, so wavelet theta
%                                        has no reference)
%   C1  battery_phasediff.csv has exactly one row per subject x condition x
%       band x decomposition x reference x pair of battery_persubject.csv,
%       with resultant_len in [0, 1]
% Tolerance: 1e-6 per subject value.
%
% Reference files are read from cfg.deriv_dir only. For each one the full
% path and last-modified date are printed; if a file of the same name also
% exists in <repo>/results/, the two are compared and FAIL is printed if
% their contents differ (CSV line endings normalised).

clear; clc;
cfg = config();
repo = fileparts(mfilename('fullpath'));
tol_val = 1e-6;
nfail = 0;

fprintf('Reading battery output from %s ...\n', cfg.deriv_dir);
B = local_read(fullfile(cfg.deriv_dir,'battery_persubject.csv'), ...
    {'subject','condition','band','measure','decomposition','reference','pair'}, {'value'});
P = local_read(fullfile(cfg.deriv_dir,'battery_phasediff.csv'), ...
    {'subject','condition','band','decomposition','reference','pair'}, {'mean_angle_deg','resultant_len'});

Rh = local_roi_avg(B, 'hilbert', 'average');
Rw = local_roi_avg(B, 'wavelet', 'average');

% ---- T1 ----
[Ref, nfail] = local_ref(cfg, repo, 'conn_results.csv', nfail);
nfail = nfail + local_cmp_roi(Ref, Rh, {'theta','beta'}, tol_val, ...
    'T1 Hilbert/average ROI-average vs conn_results.csv');

% ---- T2 ----
[Ref, nfail] = local_ref(cfg, repo, 'conn_results_betasplit.csv', nfail);
nfail = nfail + local_cmp_roi(Ref, Rh, {'beta_low','beta_high'}, tol_val, ...
    'T2 Hilbert/average ROI-average vs conn_results_betasplit.csv');

% ---- T3 ----
[Ref, nfail] = local_ref(cfg, repo, 'conn_results_wavelet.csv', nfail);
if ~isempty(Ref)
    wb = intersect(unique(Ref.band), fieldnames(cfg.battery.bands), 'stable');
    fprintf('  wavelet bands with a reference: %s\n', strjoin(wb, ', '));
    nfail = nfail + local_cmp_roi(Ref, Rw, wb, tol_val, ...
        'T3 Wavelet/average ROI-average vs conn_results_wavelet.csv');
end

% ---- C1 ----
nfail = nfail + local_check_phase(B, P);

if nfail == 0
    fprintf('\nOVERALL: PASS\n');
else
    fprintf('\nOVERALL: FAIL (%d check(s) failed)\n', nfail);
end

% ======================= local functions =======================
function T = local_read(f, textVars, numVars)
    if ~isfile(f), error('missing %s -- run s18_battery.m first', f); end
    opts = detectImportOptions(f, 'Delimiter', ',');
    opts = setvartype(opts, textVars, 'categorical');
    opts = setvartype(opts, numVars, 'double');
    T = readtable(f, opts);
end

function R = local_roi_avg(B, dec, rf)
    % ROI average over pairs, mean(...,'omitnan') as s03/s09/s16
    S = B(B.decomposition == dec & B.reference == rf, :);
    [G, subject, condition, band, measure] = findgroups(S.subject, S.condition, S.band, S.measure);
    value = splitapply(@(v) mean(v, 'omitnan'), S.value, G);
    R = table(subject, condition, band, measure, value);
end

function [T, nfail] = local_ref(cfg, repo, name, nfail)
    T = [];
    f = fullfile(cfg.deriv_dir, name);
    fprintf('\nReference: %s\n', f);
    if ~isfile(f)
        fprintf('FAIL  reference file not found\n'); nfail = nfail + 1; return;
    end
    d = dir(f);
    fprintf('  last modified: %s\n', d.date);
    alt = fullfile(repo, 'results', name);
    if isfile(alt)
        if local_same_content(f, alt)
            fprintf('PASS  identical to %s\n', alt);
        else
            fprintf('FAIL  contents differ from %s\n', alt); nfail = nfail + 1;
        end
    end
    T = readtable(f, 'TextType', 'char');
end

function same = local_same_content(a, b)
    [~,~,ext] = fileparts(a);
    if strcmpi(ext, '.mat')
        same = isequaln(load(a), load(b));
    else
        crlf = sprintf('\r\n');
        same = strcmp(strrep(fileread(a), crlf, newline), strrep(fileread(b), crlf, newline));
    end
end

function nf = local_cmp_roi(Ref, Bat, bands, tol, label)
    % Ref: subject, load, band, dwpli, plv (one row per subject x load x band)
    nf = 0;
    if isempty(Ref), fprintf('FAIL  %s: no reference loaded\n', label); nf = 1; return; end
    Ref = Ref(ismember(Ref.band, bands), :);
    kr = string(Ref.subject) + "|" + string(Ref.load) + "|" + string(Ref.band);
    meas = {'dwpli','dwPLI'; 'plv','PLV'};
    for m = 1:size(meas,1)
        Bm = Bat(Bat.measure == meas{m,2} & ismember(string(Bat.band), bands), :);
        kb = string(Bm.subject) + "|" + string(Bm.condition) + "|" + string(Bm.band);
        [tf, loc] = ismember(kr, kb);
        rv = Ref.(meas{m,1})(tf); bv = Bm.value(loc(tf));
        d = abs(rv - bv); d(isnan(rv) & isnan(bv)) = 0;
        n_missing = sum(~tf); n_extra = sum(~ismember(kb, kr));
        maxd = max([d; 0]);
        ok = n_missing == 0 && n_extra == 0 && all(d <= tol);
        fprintf('%s  %s, %s: %d rows compared, %d missing in battery, %d extra in battery, max |diff| = %.3g\n', ...
            local_pf(ok), label, meas{m,2}, sum(tf), n_missing, n_extra, maxd);
        nf = nf + ~ok;
    end
end

function nf = local_check_phase(B, P)
    kB = unique(string(B.subject) + "|" + string(B.condition) + "|" + string(B.band) + "|" + ...
        string(B.decomposition) + "|" + string(B.reference) + "|" + string(B.pair));
    kP = string(P.subject) + "|" + string(P.condition) + "|" + string(P.band) + "|" + ...
        string(P.decomposition) + "|" + string(P.reference) + "|" + string(P.pair);
    n_dup = numel(kP) - numel(unique(kP));
    n_diff = numel(setxor(kB, kP));
    r = P.resultant_len;
    n_bad_r = sum(~isnan(r) & (r < 0 | r > 1 + 1e-12));
    ok = n_dup == 0 && n_diff == 0 && n_bad_r == 0;
    fprintf('\n%s  C1 battery_phasediff.csv vs battery_persubject.csv cells: %d rows, %d duplicate, %d unmatched, %d resultant_len outside [0,1]\n', ...
        local_pf(ok), numel(kP), n_dup, n_diff, n_bad_r);
    nf = ~ok;
end

function s = local_pf(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
