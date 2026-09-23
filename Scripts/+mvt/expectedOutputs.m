function outputs = expectedOutputs(stage, day, opts)
% MVT.EXPECTEDOUTPUTS  Files a stage is expected to produce for one day.
%
% Purpose
%   One declaration of "what does this stage produce", shared by the staleness
%   check, the Makefile stamps, and the output-expectation tests. Data files
%   are named exactly; figures are given as globs, because their names encode
%   plotting parameters (direction, lane, time window) that are tunable
%   constants inside the plotting scripts.
%
% Inputs
%   stage  char, one of:
%            'gps'      assemble_data_GPS
%            'slim'     generate_data_mvt_slim
%            'full'     generate_data_mvt_full
%            'lanes'    generate_orig_dist_lanes (per-segment sidecars)
%            'lc'       extract_lane_changes_v_dist_to_av
%            'relspeed' relative_speed_histogram (pooled samples)
%            'relspeedplot' plot_relative_speed (+ binned_relative_speed)
%            'lcplot'   plotting_LC_analysis
%            'samples'  generate_data_samples
%            'fields'   generate_macroscopic_fields
%            'macro'    plot_macroscopic_fields
%            'micro'    plot_microscopic_trajectories
%            'av'       plot_AV_analysis (cross-day; day selects the folder
%                       the current code writes into)
%   day    16, 17, or 18
%   opts   (optional) options struct, passed to mvt.manifest
%
% Outputs
%   outputs  cellstr of absolute paths, possibly containing '*' globs
%
% Notes
%   'slim'/'full'/'lanes' consult mvt.manifest, so the caller learns all 24
%   expected filenames without decoding any raw data.
%
% Dependencies
%   mvt.paths, mvt.dayDir, mvt.manifest, mvt.options
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 3 || isempty(opts)
    opts = mvt.options();
end
mvt.assertDay(day);

figuresDir = mvt.dayDir('figures', day);
% The derived .mat intermediates live apart from the generated figures, so a
% results-only download can take the inputs without also taking the outputs.
analysisDir = mvt.dayDir('analysis', day);
dateTag = sprintf('202211%d', day);   % yyyyMMdd as used in figure names

switch lower(stage)
    case 'gps'
        outputs = {fullfile(mvt.dayDir('gps', day), ...
            sprintf('CIRCLES_GPS_10Hz_2022-11-%d.json', day))};

    case {'slim', 'full'}
        outDir = mvt.dayDir(lower(stage), day);
        outputs = segmentOutputs(outDir, day, opts, '');

    case 'lanes'
        % One sidecar per raw segment, named after the slim file it pairs with
        % (I-24MOTION_<date>_<time>_orig_dist_lane.mat). Derived from the
        % manifest, so the 24 names are known without decoding anything.
        %
        % These live under analysis/ rather than beside the slim JSON they
        % describe: slim/ is the released data set, and a per-segment
        % intermediate of a downstream analysis does not belong in it.
        outputs = segmentOutputs(analysisDir, day, opts, 'sidecar');

    case 'lc'
        outputs = {fullfile(analysisDir, sprintf('LC_data_%d.mat', day))};

    case 'relspeed'
        % The pooled relative-speed samples: the expensive half, so that the
        % figures can be remade without the slim tree.
        outputs = {fullfile(analysisDir, sprintf('relspeed_data_%d.mat', day))};

    case 'relspeedplot'
        % Order matters: plot_relative_speed saves against this list by index.
        outputs = { ...
            fullfile(figuresDir, sprintf('fig_relspeed_histogram_%s.pdf', dateTag)), ...
            fullfile(figuresDir, sprintf('fig_relspeed_behind_av_%s.pdf', dateTag))};

    case 'lcplot'
        % Order matters: plotting_LC_analysis saves its four figures against
        % this list by index.
        outputs = { ...
            fullfile(figuresDir, sprintf('fig_lc_exposure_%s.png', dateTag)), ...
            fullfile(figuresDir, sprintf('fig_lc_rate_merge_out_%s.png', dateTag)), ...
            fullfile(figuresDir, sprintf('fig_lc_cumulative_excess_%s.png', dateTag)), ...
            fullfile(figuresDir, sprintf('fig_lc_cumulative_excess_combined_%s.png', dateTag))};

    case 'samples'
        outputs = {fullfile(analysisDir, ...
            sprintf('samples_for_distance_analysis_%d.mat', day))};

    case 'fields'
        outputs = {fullfile(analysisDir, ...
            sprintf('fields_motion_2022-11-%d.mat', day))};

    case 'macro'
        % fig_field_<yyyyMMdd>_<direction>_<lane>_motion_<field>[_av]_nature_large.png
        outputs = {fullfile(figuresDir, ...
            sprintf('fig_field_%s_*_motion_*_nature_large.png', dateTag))};

    case 'micro'
        % fig_motion_trajectories_<yyyyMMdd>_<direction>_<lane>_lowres.png
        outputs = {fullfile(figuresDir, ...
            sprintf('fig_motion_trajectories_%s_*_lowres.png', dateTag))};

    case 'av'
        % Cross-day figures: they live in the shared figures folder, not in
        % any single day's folder, so `day` only satisfies the signature.
        p = mvt.paths();
        figuresRoot = fullfile(p.resultsDir, 'figures');
        outputs = { ...
            fullfile(figuresRoot, 'fig_2_fuel_results_*.png'), ...
            fullfile(figuresRoot, 'fig_3_fuel_results_*.png'), ...
            fullfile(figuresRoot, 'fig_SM2_vehicle_samples_counts_*.png')};

    otherwise
        error('mvt:expectedOutputs:unknownStage', ...
            ['Unknown stage ''%s''. Expected one of: gps, slim, full, ', ...
            'lanes, lc, lcplot, relspeed, relspeedplot, samples, fields, macro, ', ...
            'micro, av.'], stage);
end
end

% ---------------------------------------------------------------------------
function outputs = segmentOutputs(outDir, day, opts, kind)
% The 24 per-segment filenames, from the raw manifest when the raw data is
% present and from the product folder itself when it is not.
%
% Someone who downloaded the released trajectories without the raw inputs has
% no manifest and cannot build one: it is derived from the raw segments. Asking
% the folder what it holds lets such a tree still be described, accepted and
% reasoned about. The trade-off is that a *missing* segment cannot be noticed
% this way - but it could not be rebuilt either, so there is nothing to report.
try
    segments = mvt.manifest(day, opts);
    outputs = cell(1, numel(segments));
    for iSeg = 1:numel(segments)
        name = segments(iSeg).outputName;
        if strcmp(kind, 'sidecar')
            name = mvt.laneSidecarName(name);
        end
        outputs{iSeg} = fullfile(outDir, name);
    end
    return
catch err
    if ~strcmp(err.identifier, 'mvt:manifest:noRawFolder') && ...
            ~contains(err.message, 'Raw MOTION folder does not exist')
        rethrow(err)
    end
end

if strcmp(kind, 'sidecar')
    pattern = 'I-24MOTION_*_orig_dist_lane.mat';
else
    pattern = 'I-24MOTION_*.json';
end
listing = dir(fullfile(outDir, pattern));
listing = listing(~startsWith({listing.name}, '.'));
if isempty(listing)
    error('mvt:expectedOutputs:noManifestNoProduct', ...
        ['Cannot list the expected segments for 2022-11-%d: the raw data is ', ...
         'absent, so no manifest can be built, and %s holds no matching ', ...
         'files to describe instead.'], day, outDir);
end
outputs = cell(1, numel(listing));
for iFile = 1:numel(listing)
    outputs{iFile} = fullfile(outDir, listing(iFile).name);
end
end
