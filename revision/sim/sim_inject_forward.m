function F = sim_inject_forward(cfg, elec, labels, cache_key, force)
% sim_inject_forward -- template BEM leadfield for a REAL recording's montage.
%
%   F = sim_inject_forward(cfg, ft_raw.elec, ft_raw.label, 'ds004117')
%
% The elec struct comes from eeglab2fieldtrip applied to the subject's own
% _clean.set, so the electrode geometry is the recording's own and FieldTrip,
% not this code, does the EEGLAB-to-FieldTrip coordinate conversion. Nothing
% about the montage is assumed or supplied from outside.
%
% Returns
%   F.labels    1 x nchan, in the recording's channel order
%   F.elec      the elec struct as used (mm)
%   F.gridpos   ngrid x 3, mm
%   F.lf        nchan x (3*ngrid) free-orientation leadfield
%   F.sites     struct with .frontal and .temporal: the grid point nearest the
%               inward radial projection of cfg.site_frontal / cfg.site_temporal,
%               with a radial orientation and its nchan x 1 sensor pattern
%
% Cached per montage: the leadfield depends only on the channel set and
% positions, so it is computed once for the whole dataset.

if nargin < 5, force = false; end
if ~isfolder(cfg.cache_dir), mkdir(cfg.cache_dir); end
labels = reshape(cellstr(labels), 1, []);

cache_file = fullfile(cfg.cache_dir, sprintf('sim_inject_forward_%s.mat', cache_key));
if isfile(cache_file) && ~force
    S = load(cache_file, 'F');
    if isequal(S.F.labels, labels)
        F = S.F;
        fprintf('injection forward model: loaded cache %s (%d ch, %d grid points)\n', ...
            cache_file, numel(F.labels), size(F.gridpos,1));
        return
    end
    warning('forward cache %s was built for a different montage; recomputing.', cache_file);
end

if ~isfile(cfg.fwd.headmodel), error('headmodel not found: %s', cfg.fwd.headmodel); end
H = load(cfg.fwd.headmodel);
if ~isfield(H,'vol'), error('%s does not contain the variable "vol".', cfg.fwd.headmodel); end
headmodel = ft_convert_units(H.vol, 'mm');
elec = ft_convert_units(elec, 'mm');
if ~isfield(elec,'elecpos') || isempty(elec.elecpos)
    if isfield(elec,'chanpos') && ~isempty(elec.chanpos)
        elec.elecpos = elec.chanpos;      % some FieldTrip versions fill only chanpos
    else
        error('the recording''s elec struct carries no electrode positions.');
    end
end

elec_labels = reshape(cellstr(elec.label), 1, []);
if ~isequal(elec_labels, labels)
    [tf, loc] = ismember(upper(labels), upper(elec_labels));
    if ~all(tf)
        error('%d channel(s) of the recording have no electrode position; positions are not invented.', ...
            sum(~tf));
    end
    e = struct('label', {reshape(labels,[],1)}, 'elecpos', elec.elecpos(loc,:), ...
        'chanpos', elec.elecpos(loc,:), 'unit', 'mm');
    elec = ft_datatype_sens(e);
end

noPos = any(~isfinite(elec.elecpos), 2);
if any(noPos)
    error('%d channel(s) lack finite 3-D positions (%s). CSD and the leadfield need them.', ...
        sum(noPos), strjoin(labels(noPos), ', '));
end

scfg = [];
scfg.headmodel  = headmodel;
scfg.elec       = elec;
scfg.resolution = cfg.fwd.grid_resolution;
scfg.unit       = cfg.fwd.grid_unit;
sourcemodel = ft_prepare_sourcemodel(scfg);

lcfg = [];
lcfg.headmodel   = headmodel;
lcfg.elec        = elec;
lcfg.sourcemodel = sourcemodel;
lcfg.reducerank  = 'no';
lfs = ft_prepare_leadfield(lcfg);

inside = find(lfs.inside(:)');
ngrid = numel(inside);
nchan = numel(labels);
if ngrid < 100
    error('only %d grid points inside the brain compartment -- check the geometry.', ngrid);
end
LF = zeros(nchan, 3*ngrid);
for k = 1:ngrid
    g = lfs.leadfield{inside(k)};
    if isempty(g) || ~isequal(size(g), [nchan 3])
        error('leadfield at inside grid point %d has size [%s], expected [%d 3].', ...
            k, num2str(size(g)), nchan);
    end
    LF(:, 3*(k-1)+(1:3)) = g;
end
if ~all(isfinite(LF(:))), error('leadfield contains non-finite entries.'); end

F = struct();
F.labels    = labels;
F.elec      = elec;
F.headmodel = headmodel;
F.gridpos   = lfs.pos(inside,:);
F.lf        = LF;

centre = mean(F.gridpos, 1);
F.sites.frontal  = local_site(F, cfg.site_frontal,  cfg.inward_frac, centre);
F.sites.temporal = local_site(F, cfg.site_temporal, cfg.inward_frac, centre);
F.built_on = datestr(now, 'yyyy-mm-dd HH:MM:SS');

save(cache_file, 'F', '-v7.3');
fprintf('injection forward model: %d channels, %d inside grid points -> %s\n', ...
    nchan, ngrid, cache_file);
fprintf('  sites: %s -> grid %d (%.0f mm from the electrode), %s -> grid %d (%.0f mm)\n', ...
    cfg.site_frontal, F.sites.frontal.idx, F.sites.frontal.dist, ...
    cfg.site_temporal, F.sites.temporal.idx, F.sites.temporal.dist);
end

% ======================= local functions =======================
function S = local_site(F, label, inward_frac, centre)
    % Grid point nearest the inward radial projection of the named electrode,
    % radial orientation. Same explicit rule the scenario simulation uses.
    j = find(strcmpi(reshape(cellstr(F.elec.label),1,[]), label), 1);
    if isempty(j), error('injection site electrode "%s" is not in the montage.', label); end
    p = F.elec.elecpos(j,:);
    target = centre + inward_frac * (p - centre);
    d = vecnorm(F.gridpos - target, 2, 2);
    [~, k] = min(d);
    S = struct();
    S.label = label;
    S.idx   = k;
    S.pos   = F.gridpos(k,:);
    S.dist  = norm(S.pos - p);
    S.ori   = (S.pos - centre) / max(norm(S.pos - centre), realmin);
    S.gain  = F.lf(:, 3*(k-1)+(1:3)) * S.ori(:);
end
