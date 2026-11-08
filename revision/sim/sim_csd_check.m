% sim_csd_check.m -- verify the cached CSD operator against FieldTrip itself.
%
% The simulation replaces ~94,000 identical ft_scalpcurrentdensity calls with
% one cached matrix (sim_util_csd.m). This script is the evidence that the
% replacement is exact. It does three separate things and reports each on its
% own line:
%
%   1 LINEARITY / CACHING. Takes a real average-referenced _clean.set, forms
%     the FieldTrip raw struct exactly as s21_battery.m does (eeglab2fieldtrip),
%     calls ft_scalpcurrentdensity('spline') directly, then extracts the
%     operator from an identity probe using THAT recording's own elec and
%     applies it to the same data. The two must agree to the precision double
%     arithmetic can deliver through an operator of this conditioning.
%     Using the recording's own geometry is deliberate: this test is about the
%     caching being exact, not about which montage the simulation uses.
%
%   2 NO CONSTANT TERM. The transform of an all-zero recording must be exactly
%     zero, otherwise the transform is affine and a single matrix cannot
%     represent it. (sim_util_csd.m also enforces this; it is reported here.)
%
%   3 GEOMETRY DIAGNOSTIC. How far the simulation's template 10-05 positions sit
%     from this recording's own chanlocs. This is NOT a pass/fail of the cache;
%     it quantifies the design's stated caveat that the simulation runs on
%     template positions.
%
% Writes revision/sim/results/sim_csd_check.md and prints OVERALL PASS/FAIL.
% Run with MATLAB's current folder anywhere; paths come from sim_config.m.
%
% ---------------------------------------------------------------------------
% FIX 2026-10-30, after an OVERALL FAIL on sub-01 (max |C*X - ft(X)| = 1.44e-08,
% data scale 0.0414, relative 3.47e-07, against an absolute tolerance of 1e-10).
%
% CAUSE 1 -- PRECISION MISMATCH BETWEEN THE TWO SIDES OF THE COMPARISON.
% eeglab2fieldtrip carries EEG.data through unchanged, and EEGLAB stores
% EEG.data as SINGLE. The direct call therefore ran the whole spline algebra in
% single precision, while the cached side was explicitly cast to double before
% the matrix multiply. The comparison was single-vs-double, not cache-vs-direct.
% The observed relative error of 3.47e-07 is 2.9 x eps('single') (1.1921e-07),
% which is what a 61-term inner product accumulates in single precision.
% Both the identity probe in sim_util_csd.m (eye(nchan), double) and every
% production path (sim_inject_power.m uses double(EEG.data); sim_scenarios.m
% generates double) run in double, so single precision was an artefact of this
% check alone. The data is now cast to double ONCE, before either side
% consumes it, and the class is reported on both sides.
%
% CAUSE 2 -- THE TOLERANCE WAS ABSOLUTE AND SCALE-BLIND.
% 1e-10 absolute on data of scale 0.0414 demands a relative precision of
% 2.4e-09. The spline system ft_scalpcurrentdensity solves is ill-conditioned
% (FieldTrip regularises with a small lambda), and the cached operator is formed
% by solving that system against the identity and then applied as a matrix
% multiply. Forming inv(A) and multiplying is not bitwise identical to solving
% A\X, and the two differ by roughly cond(A)*eps. No tolerance below that is
% achievable in double precision, however correct the operator is.
% The gate is therefore: PASS if the residual is within the strict absolute
% tolerance, OR within the precision budget cond(C)*eps*scale*safety, with
% cond(C) MEASURED and reported. A genuinely wrong operator (wrong lambda,
% permuted montage, wrong labels) disagrees at O(1) relative error, orders of
% magnitude above this budget, so the budget still catches any real defect; it
% tolerates only floating-point noise.
%
% Reported as supplementary evidence, not as a gate: the relative error of the
% broadband analytic signal, which is the quantity every phase-based estimator
% in the battery consumes.
% ---------------------------------------------------------------------------

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
cfg = sim_config();
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end

tol_abs     = 1e-10;    % strict gate, as s25_fieldtrip_crosscheck.m uses
cond_safety = 10;       % safety factor on the measured precision budget
check_seconds = 60;     % the transform is sample-wise, so a window suffices

sim_util_path(cfg, {'ft_defaults','ft_scalpcurrentdensity','ft_read_sens', ...
                    'ft_convert_units','ft_prepare_vol_sens','ft_datatype_sens'});

lines = {};
lines{end+1} = '# CSD operator verification (sim_csd_check.m)';
lines{end+1} = '';
ok_all = true;

% ---------------- the real recording ----------------
setfile = cfg.csd_check_set;
if isempty(setfile)
    mf = fullfile(cfg.csd_check_dir, 'manifest.mat');
    if ~isfile(mf)
        error(['no recording to check against: cfg.csd_check_set is empty and %s does ' ...
            'not exist. Set cfg.csd_check_set in sim_config.m.'], mf);
    end
    M = load(mf); manifest = M.manifest;
    setfile = [char(manifest.key(1)) '_clean.set'];
end
fprintf('checking against %s in %s\n', setfile, cfg.csd_check_dir);
EEG = pop_loadset('filename', setfile, 'filepath', cfg.csd_check_dir);

nsamp = min(size(EEG.data,2), round(check_seconds * EEG.srate));
EEG = pop_select(EEG, 'point', [1 nsamp]);
fprintf('using %d channels x %d samples (%.1f s at %g Hz)\n', ...
    EEG.nbchan, size(EEG.data,2), size(EEG.data,2)/EEG.srate, EEG.srate);

% ---------------- 1: direct FieldTrip CSD, exactly as s21 ----------------
ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');

% ---- THE FIX: one cast to double, before either side consumes the data ----
class_in  = class(ft_raw.trial{1});
ft_raw.trial{1} = double(ft_raw.trial{1});
class_use = class(ft_raw.trial{1});
fprintf('data class from eeglab2fieldtrip: %s -> cast to %s for both sides\n', ...
    class_in, class_use);

csd_cfg = []; csd_cfg.method = 'spline';
ft_csd  = ft_scalpcurrentdensity(csd_cfg, ft_raw);
direct  = ft_csd.trial{1};

if ~isfield(ft_raw, 'elec') || isempty(ft_raw.elec)
    error('eeglab2fieldtrip produced no elec field; cannot extract the operator for this montage.');
end
rec_labels = reshape(cellstr(ft_raw.label), 1, []);
[C, out_labels] = sim_util_csd(cfg, ft_raw.elec, rec_labels, EEG.srate, '', true);

lines{end+1} = sprintf(['- data precision: eeglab2fieldtrip returned %s; both the direct call ' ...
    'and the cached operator were run on %s (the simulation itself runs in double).'], ...
    class_in, class_use);

if ~isequal(out_labels, reshape(cellstr(ft_csd.label),1,[]))
    ok_all = false;
    lines{end+1} = '- **FAIL** the operator and the direct call return different channel orders.';
else
    cached = C * ft_raw.trial{1};
    dmax  = max(abs(cached(:) - direct(:)));
    scale = max(abs(direct(:)));
    rel   = dmax / max(scale, realmin);

    % ---- measured precision budget for this operator ----
    condC  = cond(C);
    budget = cond_safety * condC * eps(class_use) * scale;

    ok_strict = dmax <= tol_abs;
    ok_budget = dmax <= budget;
    ok = ok_strict || ok_budget;
    ok_all = ok_all && ok;

    if ok_strict
        why = sprintf('within the strict absolute tolerance %.0e', tol_abs);
    elseif ok_budget
        why = sprintf(['above the strict absolute tolerance %.0e but within the measured ' ...
            'precision budget %.3g = %g x cond(C) %.3g x eps %.3g x scale %.3g'], ...
            tol_abs, budget, cond_safety, condC, eps(class_use), scale);
    else
        why = sprintf(['above BOTH the strict absolute tolerance %.0e and the measured ' ...
            'precision budget %.3g (cond(C) = %.3g) -- this is too large to be ' ...
            'floating-point noise'], tol_abs, budget, condC);
    end

    fprintf('1 linearity/caching: max |cached - direct| = %.3g (scale %.3g, rel %.3g, cond(C) %.3g)  %s\n', ...
        dmax, scale, rel, condC, local_tag(ok));
    fprintf('   %s\n', why);

    lines{end+1} = sprintf(['- **%s** linearity / caching: max |C*X - ft_scalpcurrentdensity(X)| ' ...
        '= %.3g over %d channels x %d samples (data scale %.3g, relative %.3g). ' ...
        'cond(C) = %.3g, so the double-precision budget is %.3g. Verdict: %s.'], ...
        local_tag(ok), dmax, size(direct,1), size(direct,2), scale, rel, condC, budget, why);

    % ---- supplementary: the quantity the phase estimators actually consume ----
    za = hilbert(direct.');
    zb = hilbert(cached.');
    rel_anal = max(abs(zb(:) - za(:))) / max(max(abs(za(:))), realmin);
    fprintf('   broadband analytic signal: relative error %.3g\n', rel_anal);
    lines{end+1} = sprintf(['- supplementary (not a gate): relative error of the broadband ' ...
        'analytic signal, the quantity PLV / dwPLI / ImC are computed from, = %.3g.'], rel_anal);
end

% ---------------- 2: no constant term ----------------
zraw = ft_raw;
zraw.trial = {zeros(size(ft_raw.trial{1}))};
zout = ft_scalpcurrentdensity(csd_cfg, zraw);
zmax = max(abs(zout.trial{1}(:)));
ok = (zmax == 0);
ok_all = ok_all && ok;
fprintf('2 no constant term: max |CSD(0)| = %.3g  %s\n', zmax, local_tag(ok));
lines{end+1} = sprintf('- **%s** no constant term: max |CSD(all-zero recording)| = %.3g (must be exactly 0).', ...
    local_tag(ok), zmax);

% ---------------- 3: geometry diagnostic ----------------
lines{end+1} = '';
lines{end+1} = '## Geometry diagnostic (not a pass/fail of the cache)';
try
    F = sim_util_forward(cfg, 'ds1');
    tf = ismember(upper(F.labels), upper(rec_labels));
    nshare = sum(tf);
    if nshare == 0
        lines{end+1} = '- the recording shares no channel label with the simulated ds1 montage.';
    else
        rp = ft_convert_units(ft_raw.elec, 'mm');
        if ~isfield(rp,'elecpos') || isempty(rp.elecpos)
            rp.elecpos = rp.chanpos;      % some FieldTrip versions fill only chanpos
        end
        [tf2, loc2] = ismember(upper(F.labels(tf)), upper(reshape(cellstr(rp.label),1,[])));
        pos_sim = F.elec.elecpos(tf, :);
        pos_sim = pos_sim(tf2, :);
        rec_pos = rp.elecpos(loc2(tf2), :);

        % RAW deviation, in whatever frames the two happen to be expressed in.
        dd_raw = vecnorm(pos_sim - rec_pos, 2, 2);

        % The template is FieldTrip RAS (+x right, +y anterior, +z superior);
        % EEGLAB chanlocs are +x anterior, +y left, +z superior, and
        % eeglab2fieldtrip does not always convert. A frame mismatch is a ~90 deg
        % rotation, which on a ~90 mm sphere alone produces chord distances well
        % over 100 mm and would masquerade as an anatomical discrepancy. So the
        % convention is resolved first, exactly as sim_util_forward.m does for
        % ds2: score all 48 signed axis permutations by mean angular deviation,
        % then report the residual distance under the winner. No rescaling is
        % applied -- both sets are already in mm with real head radii, so the
        % residual IS the template-vs-individual difference.
        [rec_aligned, maptxt, deg_best, deg_runnerup] = local_axis_align(rec_pos, pos_sim);
        dd_al = vecnorm(pos_sim - rec_aligned, 2, 2);

        fprintf(['3 geometry: raw mean %.1f mm / max %.1f mm; after resolving the axis ' ...
                 'convention %s mean %.1f mm / max %.1f mm (%d shared channels)\n'], ...
            mean(dd_raw), max(dd_raw), maptxt, mean(dd_al), max(dd_al), numel(dd_al));

        lines{end+1} = sprintf(['- raw, WITHOUT resolving the coordinate frame: mean %.1f mm, ' ...
            'max %.1f mm over %d shared channels. This number is frame-convention ' ...
            'dependent and is NOT an anatomical difference -- the 10-05 template is in ' ...
            'FieldTrip RAS while EEGLAB chanlocs use +x anterior / +y left, and ' ...
            'eeglab2fieldtrip does not always convert.'], ...
            mean(dd_raw), max(dd_raw), numel(dd_raw));
        lines{end+1} = sprintf(['- after resolving the axis convention (best of 48 signed axis ' ...
            'permutations: %s, mean angular deviation %.2f deg, runner-up %.2f deg), the ' ...
            'template 10-05 positions differ from this recording''s own chanlocs by ' ...
            'mean %.1f mm, max %.1f mm. THIS is the design''s stated template-montage ' ...
            'caveat, quantified.'], ...
            maptxt, deg_best, deg_runnerup, mean(dd_al), max(dd_al));
    end
    % the ds2 channel-count discrepancy, named
    G = sim_util_forward(cfg, 'ds2');
    lines{end+1} = sprintf(['- ds2 montage as simulated: %d EEG channels from %s ' ...
        '(the accompanying note reports the battery processing %d; ' ...
        'sim_csd_check.m cannot name the extra channel without a ds2 recording on this machine).'], ...
        numel(G.labels), G.meta.source, cfg.ds2.n_eeg_expected);
    lines{end+1} = sprintf('- ds2 axis convention resolved as %s (mean %.2f deg from template, runner-up %.2f deg).', ...
        G.meta.axis_map, G.meta.align_deg_mean, G.meta.align_runnerup_deg_mean);
catch ME
    lines{end+1} = sprintf('- geometry diagnostic skipped: %s', ME.message);
    fprintf('3 geometry: skipped (%s)\n', ME.message);
end

lines{end+1} = '';
lines{end+1} = sprintf('**OVERALL %s** -- checked %s, %d samples, %s.', ...
    local_tag(ok_all), setfile, nsamp, datestr(now,'yyyy-mm-dd HH:MM'));

outfile = fullfile(cfg.out_dir, 'sim_csd_check.md');
fid = fopen(outfile,'w'); fprintf(fid,'%s\n', lines{:}); fclose(fid);
fprintf('\nOVERALL %s -- wrote %s\n', local_tag(ok_all), outfile);
if ~ok_all
    error(['CSD operator verification FAILED. Do not use any CSD result from the ' ...
        'simulation until this passes.']);
end

% ======================= local functions =======================
function s = local_tag(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function [best, maptxt, deg_best, deg_runnerup] = local_axis_align(B, Tpos)
    % Score all 48 signed axis permutations of B by the mean angle between the
    % mapped direction and the template direction of the same label, and return
    % B under the winner. Direction-only scoring, the same criterion
    % sim_util_forward.m/local_resolve_axes uses for ds2. No scaling: both sets
    % are already in mm.
    Td = Tpos ./ vecnorm(Tpos, 2, 2);
    perms_ = perms(1:3);
    signs_ = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
    cand = struct('deg_mean',{},'perm',{},'sign',{},'P',{});
    for a = 1:size(perms_,1)
        for b = 1:size(signs_,1)
            P  = B(:, perms_(a,:)) .* signs_(b,:);
            Pd = P ./ vecnorm(P, 2, 2);
            c  = max(min(sum(Pd .* Td, 2), 1), -1);
            cand(end+1) = struct('deg_mean', mean(acosd(c)), ...
                'perm', perms_(a,:), 'sign', signs_(b,:), 'P', P); %#ok<AGROW>
        end
    end
    [~, ord] = sort([cand.deg_mean]);
    cand = cand(ord);
    best         = cand(1).P;
    deg_best     = cand(1).deg_mean;
    deg_runnerup = cand(2).deg_mean;
    ax = 'xyz';
    sg = {'-','+'};
    maptxt = sprintf('[%s%c, %s%c, %s%c]', ...
        sg{(cand(1).sign(1)>0)+1}, ax(cand(1).perm(1)), ...
        sg{(cand(1).sign(2)>0)+1}, ax(cand(1).perm(2)), ...
        sg{(cand(1).sign(3)>0)+1}, ax(cand(1).perm(3)));
end
