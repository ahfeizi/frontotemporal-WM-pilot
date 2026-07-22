% s00_inventory.m  --  RUN THIS FIRST.
% Loads one subject and prints channels, event types, and event fields so you
% can fill in the TODO_* placeholders in config.m. Nothing downstream is valid
% until the maintenance-onset event label and the load field/values are correct.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;

% Import BIDS (loads all subjects into ALLEEG; we inspect the first).
[~, ALLEEG] = pop_importbids(cfg.bids_dir, 'bidsevent','on', 'bidschanloc','on');
EEG = eeg_checkset(ALLEEG(1));

fprintf('\n========== ds004117 INVENTORY ==========\n');
fprintf('Subjects loaded: %d\n', numel(ALLEEG));
fprintf('Sampling rate  : %g Hz\n', EEG.srate);
fprintf('Duration       : %.1f s\n', EEG.pnts/EEG.srate);

fprintf('\n--- CHANNELS (%d) ---\n', EEG.nbchan);
disp({EEG.chanlocs.labels});

fprintf('\n--- UNIQUE EVENT TYPES ---\n');
types = {EEG.event.type};
u = unique(types);
for i = 1:numel(u)
    fprintf('  %-24s n=%d\n', u{i}, sum(strcmp(types, u{i})));
end

fprintf('\n--- EVENT FIELDS ---\n');
disp(fieldnames(EEG.event));

fprintf('\n--- SAMPLE EVENTS (first 15) ---\n');
struct2table(EEG.event(1:min(15,numel(EEG.event))))

fprintf(['\nACTION: set in config.m ->\n' ...
         '  cfg.ev_maint_lock = <the maintenance/retention onset label>\n' ...
         '  cfg.load_field    = <the field holding memory load>\n' ...
         '  cfg.load_low / cfg.load_high = <the actual load values>\n' ...
         'Also verify roi_frontal/roi_temporal exist in the channel list above.\n']);
