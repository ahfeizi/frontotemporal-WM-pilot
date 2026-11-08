function V = sim_inject_battery(cfg, runs, win, wlen, fr, te, do_wavelet)
% sim_inject_battery -- the battery on ONE arm of the injection experiment.
%
%   V = sim_inject_battery(cfg, runs, win, wlen, fr, te, do_wavelet)
%
% runs : struct array, one element per run of the subject, with fields
%        .srate, .locks (lock samples of the retained high-load trials),
%        .average and .csd (ROI-row continuous data for this arm)
% win, wlen : window in samples, as the battery computes them
% fr, te    : row indices of the frontal and temporal ROI rows
%
% Returns V.(ref).(decomp).(band).(measure) = 1 x np, and V.ntrials_used.
%
% The code path is s18_battery.m's: continuous band-pass + Hilbert (via
% util_bandpass.m) or unit-energy complex Morlet at 5 cycles, trials pooled
% ACROSS RUNS exactly as s03_connectivity.m pools them, the maintenance window
% cut after the continuous transform, out-of-bounds trials dropped, cells with
% fewer than cfg.min_trials segments left as NaN, and PLV/dwPLI/ImC from
% util_conn_pairs.m. s18's local functions are copied verbatim because MATLAB
% local functions cannot be shared between files -- not re-derived.

if nargin < 7, do_wavelet = false; end
bands = reshape(cfg.bands_used, 1, []);
meas  = cfg.measures;
np    = numel(fr) * numel(te);

V = struct();
V.ntrials_used = 0;

for r = 1:numel(cfg.refs)
    rf = cfg.refs{r};

    % ---------------- band-pass + Hilbert (as s18/s03) ----------------
    for b = 1:numel(bands)
        bn = bands{b}; band = cfg.bands.(bn);
        seg = local_pool(runs, rf, win, wlen, ...
            @(X, srate) local_hilbert(X, srate, band, cfg.filt_order_factor));
        for m = 1:numel(meas)
            V.(rf).hilbert.(bn).(meas{m}) = nan(1, np);
        end
        kept = local_ntr(seg);
        V.ntrials_used = kept;
        if kept < cfg.min_trials, continue; end
        [plv, dwpli, imc] = util_conn_pairs(seg, fr, te);
        V.(rf).hilbert.(bn).PLV   = plv;
        V.(rf).hilbert.(bn).dwPLI = dwpli;
        V.(rf).hilbert.(bn).ImC   = imc;
    end

    if ~do_wavelet, continue; end

    % ---------------- Morlet wavelet (as s18, ported from s16) ----------------
    all_freqs = [];
    for b = 1:numel(bands), all_freqs = [all_freqs cfg.wavelet_freqs.(bands{b})]; end %#ok<AGROW>
    all_freqs = unique(all_freqs);
    nf = numel(all_freqs);
    W = struct('ok', false(1,nf), 'plv', nan(nf,np), 'dwpli', nan(nf,np), 'imc', nan(nf,np));
    for f = 1:nf
        seg = local_pool(runs, rf, win, wlen, ...
            @(X, srate) local_morlet(X, srate, all_freqs(f), cfg.n_cycles));
        if local_ntr(seg) < cfg.min_trials, continue; end
        [W.plv(f,:), W.dwpli(f,:), W.imc(f,:)] = util_conn_pairs(seg, fr, te);
        W.ok(f) = true;
    end
    for b = 1:numel(bands)
        bn = bands{b};
        fb = find(ismember(all_freqs, cfg.wavelet_freqs.(bn)));
        for m = 1:numel(meas)
            V.(rf).wavelet.(bn).(meas{m}) = nan(1, np);
        end
        used = fb(W.ok(fb));
        nu = numel(used);
        if nu == 0, continue; end
        % per-pair average over bins, sum / n_used as s16 (no omitnan)
        V.(rf).wavelet.(bn).PLV   = sum(W.plv(used,:),1)   / nu;
        V.(rf).wavelet.(bn).dwPLI = sum(W.dwpli(used,:),1) / nu;
        V.(rf).wavelet.(bn).ImC   = sum(W.imc(used,:),1)   / nu;
    end
end
end

% ======================= local functions (copied from s18_battery.m) =======
function Z = local_hilbert(X, srate, band, order_factor)
    Xf = util_bandpass(X, srate, band, order_factor);
    Z  = hilbert(Xf.').';
end

function Z = local_morlet(data, srate, fc, n_cycles)
    sigma_t = n_cycles / (2*pi*fc);
    t = -3*sigma_t : 1/srate : 3*sigma_t;
    wavelet = exp(2i*pi*fc.*t) .* exp(-t.^2/(2*sigma_t^2));
    wavelet = wavelet / sqrt(sum(abs(wavelet).^2));   % unit energy
    Z = zeros(size(data));
    for ch = 1:size(data,1)
        Z(ch,:) = conv(data(ch,:), wavelet, 'same');
    end
end

function seg = local_pool(runs, rf, win, wlen, transform)
    % Maintenance-window segments pooled across the subject's runs, as s03/s18.
    seg = [];
    for k = 1:numel(runs)
        Z = transform(runs(k).(rf), runs(k).srate);
        seg = cat(3, seg, local_segments(Z, runs(k).locks, win, wlen));
    end
end

function seg = local_segments(Z, locks, win, wlen)
    seg = complex(zeros(size(Z,1), wlen, numel(locks))); kept = 0;
    for t = 1:numel(locks)
        a = locks(t) + win(1); z = a + wlen - 1;
        if a < 1 || z > size(Z,2), continue; end
        kept = kept + 1;
        seg(:,:,kept) = Z(:, a:z);
    end
    seg = seg(:,:,1:kept);
end

function n = local_ntr(seg)
    if isempty(seg), n = 0; else, n = size(seg,3); end
end
