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
%   'slim'/'full' consult mvt.manifest, so the caller learns all 24 expected
%   filenames without decoding any raw data.
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
dateTag = sprintf('202211%d', day);   % yyyyMMdd as used in figure names

switch lower(stage)
    case 'gps'
        outputs = {fullfile(mvt.dayDir('gps', day), ...
            sprintf('CIRCLES_GPS_10Hz_2022-11-%d.json', day))};

    case {'slim', 'full'}
        segments = mvt.manifest(day, opts);
        outDir = mvt.dayDir(lower(stage), day);
        outputs = cell(1, numel(segments));
        for iSeg = 1:numel(segments)
            outputs{iSeg} = fullfile(outDir, segments(iSeg).outputName);
        end

    case 'samples'
        outputs = {fullfile(figuresDir, ...
            sprintf('samples_for_distance_analysis_%d.mat', day))};

    case 'fields'
        outputs = {fullfile(figuresDir, ...
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
            'samples, fields, macro, micro, av.'], stage);
end
end
