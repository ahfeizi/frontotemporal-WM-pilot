% s11_accuracy.m -- Reviewer-requested check (C3): behavioral accuracy by load.
%
% For each subject (pooled across runs), joins:
%   - WM (maintenance onset) lock events -> trial number + load condition
%   - feedback events (task_role)        -> trial number + correct/incorrect/bad
% to compute per-subject, per-load accuracy (proportion correct), then runs
% the same sign-flip permutation test used elsewhere (s04-style) on the
% low-vs-high accuracy contrast.
%
% Uses ONLY the already-preprocessed _clean.set files: task_role, trial, and
% memory_cond survive filtering/ASR/ICA/re-referencing unchanged (s02_epoch_new.m
% already relies on this), so no re-import from raw BIDS is needed.
%
% Trial numbers are per-RUN, not globally unique, so lock/feedback events are
% joined using a compound key "<dataset_key>_<trial>" to avoid cross-run collisions
% when a subject's runs are pooled.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

subjects = unique(manifest.subject, 'stable');
rows = {};  % subject, load, n_correct, n_incorrect, n_bad, n_scored, accuracy

for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);

    loadMap = containers.Map('KeyType','char','ValueType','char');  % "key_trial" -> 'low'/'high'
    roleMap = containers.Map('KeyType','char','ValueType','char');  % "key_trial" -> task_role

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        ev = EEG.event;
        if ~isfield(ev,'trial') || ~isfield(ev,cfg.role_field)
            warning('%s: missing trial/%s fields, skipping', key, cfg.role_field);
            continue;
        end

        % --- lock events (WM onset) -> load per trial ---
        locks = find(strcmp({ev.type}, cfg.ev_maint_lock));
        for L = locks
            trN = ev(L).trial;
            if ~isnumeric(trN) || isempty(trN) || isnan(trN), continue; end
            lv = ''; if isfield(ev, cfg.load_field), lv = ev(L).(cfg.load_field); end
            lab = '';
            if ismember(string(lv), string(cfg.load_low)),  lab = 'low';  end
            if ismember(string(lv), string(cfg.load_high)), lab = 'high'; end
            if isempty(lab), continue; end
            loadMap(sprintf('%s_%d', key, trN)) = lab;
        end

        % --- feedback events -> role per trial ---
        roles = {ev.(cfg.role_field)};
        for e = 1:numel(ev)
            trN = ev(e).trial;
            if ~isnumeric(trN) || isempty(trN) || isnan(trN), continue; end
            r = roles{e};
            if any(strcmp(r, {cfg.role_feedback_correct, cfg.role_feedback_incorrect, cfg.role_bad}))
                roleMap(sprintf('%s_%d', key, trN)) = r;
            end
        end
    end

    % --- join on compound trial key ---
    tkeys = loadMap.keys;
    nC = struct('low',0,'high',0);
    nI = struct('low',0,'high',0);
    nB = struct('low',0,'high',0);
    for i = 1:numel(tkeys)
        tk = tkeys{i};
        lab = loadMap(tk);
        if ~isKey(roleMap, tk), continue; end   % no scored feedback for this trial
        r = roleMap(tk);
        switch r
            case cfg.role_feedback_correct,   nC.(lab) = nC.(lab) + 1;
            case cfg.role_feedback_incorrect, nI.(lab) = nI.(lab) + 1;
            case cfg.role_bad,                nB.(lab) = nB.(lab) + 1;
        end
    end

    for c = {'low','high'}
        lab = c{1};
        nCorrect = nC.(lab); nIncorrect = nI.(lab); nBad = nB.(lab);
        nScored = nCorrect + nIncorrect;   % accuracy excludes 'bad' (aborted/no-response) trials
        acc = NaN; if nScored > 0, acc = nCorrect / nScored; end
        rows(end+1,:) = {char(subj), lab, nCorrect, nIncorrect, nBad, nScored, acc}; %#ok<SAGROW>
    end
    fprintf('done %s\n', subj);
end

R = cell2table(rows, 'VariableNames', {'subject','load','n_correct','n_incorrect', ...
    'n_bad','n_scored','accuracy'});
writetable(R, fullfile(cfg.deriv_dir,'accuracy_by_subject_load.csv'));

% ---- group-level descriptives + paired permutation test (same style as s04) ----
rng(1);
subs = unique(R.subject);
lo = nan(numel(subs),1); hi = nan(numel(subs),1);
for i = 1:numel(subs)
    rl = R(strcmp(R.subject,subs{i}) & strcmp(R.load,'low'), :);
    rh = R(strcmp(R.subject,subs{i}) & strcmp(R.load,'high'), :);
    if ~isempty(rl), lo(i) = rl.accuracy(1); end
    if ~isempty(rh), hi(i) = rh.accuracy(1); end
end
ok = ~isnan(lo) & ~isnan(hi);
lo = lo(ok); hi = hi(ok); n = sum(ok);
d = hi - lo;

tobs = mean(d)/(std(d)/sqrt(n));
np_ = cfg.stat.n_perm; tnull = zeros(np_,1);
for pp = 1:np_
    dp = d .* sign(rand(n,1)-0.5);
    tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
end
pval = mean(abs(tnull) >= abs(tobs));
dz = mean(d)/std(d);

fprintf('\nAccuracy (proportion correct), N = %d\n', n);
fprintf('  Low load : mean = %.3f, SD = %.3f\n', mean(lo), std(lo));
fprintf('  High load: mean = %.3f, SD = %.3f\n', mean(hi), std(hi));
fprintf('  High-Low : dz = %.2f, t = %.2f, p = %.4f\n', dz, tobs, pval);
if mean(lo) > 0.95 && mean(hi) > 0.95
    fprintf('  NOTE: both conditions above 95%% accuracy -- possible ceiling effect.\n');
end

save(fullfile(cfg.deriv_dir,'accuracy_stats.mat'), 'R','lo','hi','dz','tobs','pval','n');
fprintf('\nSaved accuracy_by_subject_load.csv and accuracy_stats.mat\n');
