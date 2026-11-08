% F_pair_hits.m
% For every pair that passed within-cell FWER (p_fwer < 0.05) in any of the
% 48 cells, report: measure, band, decomposition, reference, t, p_fwer,
% and the t-statistic and (uncorrected) p for that same pair in the other 47 cells.
%
% Reads:  results/unified_contrasts_ds004117.csv (pair-level rows, level=='pair')
% Writes: revision/F_pair_hits.csv
%         revision/F_pair_hits.md

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();

csv_in = fullfile(fileparts(repo), 'results', 'unified_contrasts_ds004117.csv');
if ~isfile(csv_in), error('F_pair_hits:missing', 'Cannot find %s', csv_in); end

fprintf('Reading %s ...\n', csv_in);
opts = detectImportOptions(csv_in, 'Delimiter', ',');
opts = setvartype(opts, {'dataset','contrast','level','band','measure','decomposition', ...
    'reference','pair','perm_method'}, 'char');
opts = setvartype(opts, {'n','dz','dz_ci_lo','dz_ci_hi','t','df','p','p_perm', ...
    'p_fwer','bf10','bf01','p_tost'}, 'double');
T = readtable(csv_in, opts);

% Keep only pair-level rows
Tp = T(strcmp(T.level, 'pair'), :);
fprintf('Pair-level rows: %d (expected 48 x 30 = %d)\n', height(Tp), 48*30);

% Find hits: pairs passing within-cell FWER < 0.05
hit_mask = ~isnan(Tp.p_fwer) & Tp.p_fwer < 0.05;
Thit = Tp(hit_mask, :);
fprintf('Pairs passing within-cell FWER < 0.05: %d\n', height(Thit));

if height(Thit) == 0
    fprintf('No pair hits found. Writing empty table.\n');
    records = table('Size', [0 10], 'VariableTypes', repmat({'double'}, 1, 10), ...
        'VariableNames', {'hit_band','hit_measure','hit_decomp','hit_ref', ...
                          'hit_pair','hit_t','hit_p_fwer', ...
                          'other_band_measure_decomp_ref','other_t','other_p'});
    writetable(records, fullfile(repo, 'F_pair_hits.csv'));
else
    % For each hit pair, collect t and p in the other 47 cells
    all_cells = unique(Tp(:, {'band','measure','decomposition','reference'}), 'rows', 'stable');
    records = {};
    for h = 1:height(Thit)
        hit = Thit(h, :);
        pair_name = hit.pair{1};
        % Hit cell details
        row0 = {hit.band{1}, hit.measure{1}, hit.decomposition{1}, hit.reference{1}, ...
                pair_name, hit.t, hit.p_fwer, ...
                sprintf('%s__%s__%s__%s', hit.band{1}, hit.measure{1}, hit.decomposition{1}, hit.reference{1}), ...
                hit.t, hit.p};
        records{end+1} = row0; %#ok<SAGROW>
        % Other 47 cells
        for c = 1:height(all_cells)
            bn = all_cells.band{c}; me = all_cells.measure{c};
            dc = all_cells.decomposition{c}; rf = all_cells.reference{c};
            if strcmp(bn, hit.band{1}) && strcmp(me, hit.measure{1}) && ...
               strcmp(dc, hit.decomposition{1}) && strcmp(rf, hit.reference{1})
                continue; % skip the hit cell itself
            end
            other = Tp(strcmp(Tp.band, bn) & strcmp(Tp.measure, me) & ...
                       strcmp(Tp.decomposition, dc) & strcmp(Tp.reference, rf) & ...
                       strcmp(Tp.pair, pair_name), :);
            if isempty(other)
                t_o = NaN; p_o = NaN;
            else
                t_o = other.t(1); p_o = other.p(1);
            end
            cell_label = sprintf('%s__%s__%s__%s', bn, me, dc, rf);
            records{end+1} = {hit.band{1}, hit.measure{1}, hit.decomposition{1}, hit.reference{1}, ...
                              pair_name, hit.t, hit.p_fwer, cell_label, t_o, p_o}; %#ok<SAGROW>
        end
    end
    varnames = {'hit_band','hit_measure','hit_decomp','hit_ref','hit_pair', ...
                'hit_t','hit_p_fwer','other_cell','other_t','other_p'};
    out_table = cell2table(vertcat(records{:}), 'VariableNames', varnames);
    writetable(out_table, fullfile(repo, 'F_pair_hits.csv'));
    fprintf('Wrote %d rows to F_pair_hits.csv\n', height(out_table));
end

% Markdown summary
md = {};
md{end+1} = '# F. Pair-level hits table';
md{end+1} = '';
md{end+1} = sprintf('Threshold: within-cell FWER (max-t sign-flip) p < 0.05.');
md{end+1} = sprintf('Source: unified_contrasts_ds004117.csv, level==pair (%d rows).', height(Tp));
md{end+1} = sprintf('Pairs passing FWER in any cell: **%d**', height(Thit));
md{end+1} = '';
if height(Thit) > 0
    for h = 1:height(Thit)
        hit = Thit(h,:);
        md{end+1} = sprintf('- %s | %s | %s | %s | pair %s | t = %.4f | p_FWER = %.4f', ...
            hit.band{1}, hit.measure{1}, hit.decomposition{1}, hit.reference{1}, ...
            hit.pair{1}, hit.t, hit.p_fwer);
    end
else
    md{end+1} = 'No hits found.';
end
md{end+1} = '';
md{end+1} = 'Columns in F_pair_hits.csv: hit_band, hit_measure, hit_decomp, hit_ref, hit_pair,';
md{end+1} = '  hit_t, hit_p_fwer, other_cell (band__measure__decomp__ref), other_t, other_p.';
md{end+1} = 'One row for the hit cell itself (other_cell == hit cell) plus 47 rows (other cells),';
md{end+1} = 'per hit pair.';
fid = fopen(fullfile(repo, 'F_pair_hits.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Done.\n');
