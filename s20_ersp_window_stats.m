% s20_ersp_window_stats.m  --  Per-subject FCz ERSP window means + regression check.
% Recomputes the per-subject ERSPs EXACTLY as s07_grandavg_ersp.m (same channel,
% epoch, baseline, freqs, cycles, trial selection, loads collapsed) and, from each
% subject's ERSP (dB, before the grand average), takes the mean over:
%   theta 4-7 Hz,  window 0.0-0.5 s
%   alpha 8-12 Hz, window 0.5-2.5 s
% Window edges are snapped to the nearest available newtimef freqs/times; the bins
% actually used are printed.
% Output: cfg.deriv_dir/ersp_window_stats.csv  (subject, band, window, ersp_db)
% Regression check: mean over subjects of the per-subject ERSPs must equal 'grand'
% in grandavg_ersp.mat (written by s07) to 1e-6 -> prints PASS/FAIL.
% s07 and all other existing files are left untouched.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

% ---- identical to s07_grandavg_ersp.m ----
chan      = 'FCz';                 % frontal-midline
epoch_s   = [-1.0 3.5];            % around WM onset (s)
baseline  = [-500 0];             % ms, pre-dash
freqrange = [2 35];
cycles    = [3 0.5];

% ---- window definitions (s20) ----
win.name   = {'theta',  'alpha'};
win.f_hz   = {[4 7],    [8 12]};
win.t_s    = {[0.0 0.5], [0.5 2.5]};
win.label  = {'0.0-0.5s', '0.5-2.5s'};
tol        = 1e-6;

subs = unique(manifest.subject,'stable');
GA = []; times = []; freqs = []; usedN = 0; subjList = {};

for su = 1:numel(subs)
    subj = subs(su);
    keys = manifest.key(manifest.subject == subj);

    % gather single-channel maintenance epochs across this subject's runs
    data = []; eplen = []; srate = [];
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        srate = EEG.srate;
        ep = round(epoch_s * srate); if isempty(eplen), eplen = ep(2)-ep(1)+1; end
        ch = find(strcmpi({EEG.chanlocs.labels}, chan), 1);
        sel = trials([trials.good]);          % collapse load conditions
        for t = 1:numel(sel)
            a = sel(t).lock_sample + ep(1); b = a + eplen - 1;
            if a < 1 || b > EEG.pnts, continue; end
            data = [data, EEG.data(ch, a:b).']; %#ok<AGROW>  % column per trial
        end
    end
    ntr = size(data,2);
    if ntr < 5, warning('%s: only %d trials, skipping', subj, ntr); continue; end

    [ersp,~,~,t_ms,f_hz] = newtimef(data, eplen, epoch_s*1000, srate, cycles, ...
        'baseline', baseline, 'freqs', freqrange, ...
        'plotersp','off','plotitc','off','verbose','off');

    if isempty(GA)
        GA = zeros(size(ersp,1), size(ersp,2), 0); times = t_ms; freqs = f_hz;
    end
    if ~isequal(size(ersp,1),numel(freqs)) || ~isequal(size(ersp,2),numel(times))
        warning('%s: ERSP grid mismatch, skipping', subj); continue;
    end
    GA(:,:,end+1) = ersp; %#ok<AGROW>
    usedN = usedN + 1; subjList{end+1} = char(subj); %#ok<AGROW>
    fprintf('ERSP done %s  (%d trials)\n', subj, ntr);
end
% ---- end of s07 copy ----

if usedN == 0, error('No subjects with a valid ERSP.'); end

% ---- resolve window bins on the newtimef grid (nearest available) ----
nW = numel(win.name);
fIdx = cell(1,nW); tIdx = cell(1,nW);
fprintf('\nWindow bins used (nearest available on the newtimef grid):\n');
for w = 1:nW
    fIdx{w} = local_span(freqs, win.f_hz{w});
    tIdx{w} = local_span(times, win.t_s{w} * 1000);      % s -> ms (newtimef times)
    fprintf('  %-5s requested %g-%g Hz, %g-%g s\n', win.name{w}, ...
        win.f_hz{w}(1), win.f_hz{w}(2), win.t_s{w}(1), win.t_s{w}(2));
    fprintf('        freqs used (%d bins, Hz): %s\n', numel(fIdx{w}), ...
        strjoin(arrayfun(@(x) sprintf('%.3f',x), freqs(fIdx{w}), 'UniformOutput',false), ', '));
    fprintf('        times used (%d bins, ms): %.1f ... %.1f\n', numel(tIdx{w}), ...
        times(tIdx{w}(1)), times(tIdx{w}(end)));
    fprintf('        times used (ms): %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f',x), times(tIdx{w}), 'UniformOutput',false), ', '));
end

% ---- per-subject window means (dB, before grand average) ----
nRow = usedN * nW;
c_subj = cell(nRow,1); c_band = cell(nRow,1); c_win = cell(nRow,1); v_db = nan(nRow,1);
r = 0;
for s = 1:usedN
    for w = 1:nW
        r = r + 1;
        blk = GA(fIdx{w}, tIdx{w}, s);
        c_subj{r} = subjList{s};
        c_band{r} = win.name{w};
        c_win{r}  = win.label{w};
        v_db(r)   = mean(blk(:));
    end
end
R = table(c_subj, c_band, c_win, v_db, ...
    'VariableNames', {'subject','band','window','ersp_db'});
writetable(R, fullfile(cfg.deriv_dir,'ersp_window_stats.csv'));
fprintf('\nSaved ersp_window_stats.csv  (%d subjects x %d windows)\n', usedN, nW);

% ---- regression check vs s07 grandavg_ersp.mat ----
G = load(fullfile(cfg.deriv_dir,'grandavg_ersp.mat'), 'grand','times','freqs','subjList','usedN');
grand_s20 = mean(GA, 3);
pass = true;
if ~isequal(size(grand_s20), size(G.grand))
    fprintf('  grid size mismatch: s20 %s vs s07 %s\n', mat2str(size(grand_s20)), mat2str(size(G.grand)));
    pass = false; maxdiff = Inf;
else
    maxdiff = max(abs(grand_s20(:) - G.grand(:)));
    if ~(maxdiff <= tol), pass = false; end
    if max(abs(times(:) - G.times(:))) > tol || max(abs(freqs(:) - G.freqs(:))) > tol
        fprintf('  times/freqs grid differs from s07\n'); pass = false;
    end
end
if usedN ~= G.usedN || ~isequal(subjList, G.subjList)
    fprintf('  subject set differs from s07 (s20 N=%d, s07 N=%d)\n', usedN, G.usedN);
    pass = false;
end
if pass, verdict = 'PASS'; else, verdict = 'FAIL'; end
fprintf('\nRegression check (mean of per-subject ERSP == s07 grand, tol %g): %s  (max |diff| = %.3g dB)\n', ...
    tol, verdict, maxdiff);

% =====================================================================
function idx = local_span(grid, lim)
% Indices of grid between the grid points nearest to lim(1) and lim(2), inclusive.
grid = grid(:).';
[~, i1] = min(abs(grid - lim(1)));
[~, i2] = min(abs(grid - lim(2)));
idx = min(i1,i2):max(i1,i2);
end
