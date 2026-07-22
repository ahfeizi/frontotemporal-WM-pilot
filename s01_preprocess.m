% s01_preprocess.m  --  Per-DATASET preprocessing, saving CONTINUOUS cleaned data.
% ds004117 loads as ~85 datasets = ~24 subjects x several runs. We therefore
%   * save each dataset under a UNIQUE key (no overwrite of a subject's runs)
%   * record key -> subject in a manifest, so later stages aggregate to the
%     SUBJECT level (N = subjects, never N = datasets).
% Data is kept continuous; ASR runs in correction mode (no discontinuities).
%
% NOTE ON RUNTIME: ICA on every dataset can take hours. Use cfg.subset_n for a
% fast first pass, then set it to [] for the final run.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
if ~exist(cfg.deriv_dir,'dir'), mkdir(cfg.deriv_dir); end

[~, ALLEEG] = pop_importbids(cfg.bids_dir, 'bidsevent','on', 'bidschanloc','on');

nAll = numel(ALLEEG);
nDo  = nAll; if ~isempty(cfg.subset_n), nDo = min(cfg.subset_n, nAll); end
fprintf('Datasets available: %d | processing: %d\n', nAll, nDo);

manifest = table('Size',[0 3], 'VariableTypes',{'string','string','string'}, ...
                 'VariableNames',{'key','subject','clean_file'});

for s = 1:nDo
    EEG = eeg_checkset(ALLEEG(s), 'loaddata');   % load this dataset's data only
    subj = EEG.subject; if isempty(subj), subj = sprintf('s%02d', s); end
    if ~startsWith(subj,'sub-'), subj = ['sub-' subj]; end
    key = sprintf('%s_d%03d', subj, s);          % unique per dataset
    fprintf('\n=== %s  (%d/%d) ===\n', key, s, nDo);

    if EEG.srate ~= cfg.resample_hz, EEG = pop_resample(EEG, cfg.resample_hz); end
    EEG = pop_eegfiltnew(EEG, 'locutoff', cfg.hp_hz);
    EEG = pop_eegfiltnew(EEG, 'hicutoff', cfg.lp_hz);
    if ~isempty(cfg.line_hz)
        EEG = pop_cleanline(EEG, 'linefreqs', cfg.line_hz, 'bandwidth', 2);
    end

    EEG.urchanlocs = EEG.chanlocs;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5, 'ChannelCriterion',0.8, ...
        'LineNoiseCriterion',4, 'Highpass','off', 'BurstCriterion',20, ...
        'WindowCriterion','off', 'BurstRejection','off', 'Distance','Euclidian');

    dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    EEG = pop_runica(EEG, 'icatype','runica', 'extended',1, 'pca', dataRank);
    EEG = pop_iclabel(EEG, 'default');
    % ICLabel order: Brain Muscle Eye Heart LineNoise ChannelNoise Other
    EEG = pop_icflag(EEG, [NaN NaN; 0.8 1; 0.8 1; 0.8 1; 0.8 1; NaN NaN; NaN NaN]);
    EEG = pop_subcomp(EEG, [], 0);

    EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
    EEG = pop_reref(EEG, []);

    pop_saveset(EEG, 'filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    manifest = [manifest; {string(key), string(subj), string([key '_clean.set'])}]; %#ok<AGROW>

    ALLEEG(s).data = []; clear EEG;
end

writetable(manifest, fullfile(cfg.deriv_dir,'manifest.csv'));
save(fullfile(cfg.deriv_dir,'manifest.mat'),'manifest');
fprintf('\nPreprocessing complete: %d datasets across %d unique subjects.\n', ...
        height(manifest), numel(unique(manifest.subject)));
fprintf('(If unique subjects << datasets, the 85 were runs, as expected.)\n');
