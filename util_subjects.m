function subs = util_subjects(cfg)
% Return subject IDs: explicit list from config, else auto-discover sub-* dirs.
if ~isempty(cfg.subjects), subs = cfg.subjects; return; end
d = dir(fullfile(cfg.bids_dir, 'sub-*'));
subs = {d([d.isdir]).name};
end
