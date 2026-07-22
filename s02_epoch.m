% s02_epoch.m  --  Build the maintenance-window trial list for each dataset.
% Reads the manifest from s01. For each clean dataset, produces `trials`:
%   .lock_sample  sample of the 'WM' (dash / maintenance onset) event
%   .load         'low' | 'high'   (from memory_cond)
%   .good         logical, passes amplitude check AND (if cfg.correct_only)
%                 belongs to a correct-response trial
% Saved as <key>_trials.mat alongside the clean data.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

for r = 1:height(manifest)
    key = char(manifest.key(r));
    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    srate = EEG.srate;
    win = round(cfg.maint_window * srate);
    roi = util_chanidx(EEG, [cfg.roi_frontal cfg.roi_temporal]);
    ev  = EEG.event;

    % trials with a correct response (match on the 'trial' field)
    % --- Trial inclusion from task_role (uniform across all subjects) ---
    roles = {ev.task_role};
    trN = nan(1, numel(ev));
    for e = 1:numel(ev)
        v = ev(e).trial;
        if isnumeric(v) && ~isempty(v), trN(e) = v; end
    end
    badTr  = unique(trN(strcmp(roles, cfg.role_bad)));
    incTr  = unique(trN(strcmp(roles, cfg.role_feedback_incorrect)));
    corrTr = unique(trN(strcmp(roles, cfg.role_feedback_correct)));
    if cfg.correct_only
        includeTr = setdiff(corrTr, [badTr incTr]);   % correct, not bad, not wrong
    else
        includeTr = setdiff(unique(trN(~isnan(trN))), badTr);  % all but bad trials
    end

    locks = find(strcmp({ev.type}, cfg.ev_maint_lock));
    trials = struct('lock_sample',{},'load',{},'good',{});

    for L = locks
        ls = round(ev(L).latency);
        a = ls + win(1); b = ls + win(2);
        if a < 1 || b > EEG.pnts, continue; end

        lv = ''; if isfield(ev, cfg.load_field), lv = ev(L).(cfg.load_field); end
        lab = '';
        if ismember(string(lv), string(cfg.load_low)),  lab = 'low';  end
        if ismember(string(lv), string(cfg.load_high)), lab = 'high'; end
        if isempty(lab), continue; end

        good = max(abs(EEG.data(roi, a:b)), [], 'all') < cfg.amp_thresh;
        good = good && ismember(ev(L).trial, includeTr);
        trials(end+1) = struct('lock_sample',ls, 'load',lab, 'good',good); %#ok<SAGROW>
    end

    save(fullfile(cfg.deriv_dir,[key '_trials.mat']),'trials');
    fprintf('%s: %2d trials (good: low=%d high=%d)\n', key, numel(trials), ...
        sum([trials.good] & strcmp({trials.load},'low')), ...
        sum([trials.good] & strcmp({trials.load},'high')));
end
