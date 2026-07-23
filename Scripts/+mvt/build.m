function build(stage, day, varargin)
% MVT.BUILD  Run one pipeline stage for one day (the entry point make calls).
%
% Purpose
%   A single, uniform command line for every stage, so the Makefile does not
%   have to know each function's signature:
%     matlab -batch "mvt.build('slim', 17)"
%     matlab -batch "mvt.build('slim', 17, 'Shard', [2 4])"
%     matlab -batch "mvt.build('av', [], 'Days', [16 17 18])"
%
% Inputs
%   stage  'gps' | 'slim' | 'full' | 'samples' | 'fields' | 'macro' | 'micro' | 'av'
%   day    16, 17, or 18; ignored (may be []) for the cross-day 'av' stage
%   ...    options struct and/or name/value pairs (see mvt.options)
%
% Outputs
%   (none; the stage writes its own files and decides for itself whether any
%   work is needed - every stage performs its own staleness check)
%
% Dependencies
%   mvt.options, mvt.log, and the stage functions in Scripts/
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
stage = lower(char(stage));

if ~strcmp(stage, 'av')
    mvt.assertDay(day);
end

mvt.log(opts, 'stage %s%s starting', stage, dayTag(stage, day));
runTimer = tic;

switch stage
    case 'gps'
        assemble_data_GPS(day, opts);
    case 'slim'
        generate_data_mvt_slim(day, opts);
    case 'full'
        generate_data_mvt_full(day, opts);
    case 'samples'
        generate_data_samples(day, opts);
    case 'fields'
        generate_macroscopic_fields(day, opts);
    case 'macro'
        plot_macroscopic_fields(day, opts);
    case 'micro'
        plot_microscopic_trajectories(day, opts);
    case 'av'
        plot_AV_analysis(opts.Days, opts);
    otherwise
        error('mvt:build:unknownStage', ...
            ['Unknown stage ''%s''. Expected one of: gps, slim, full, ', ...
            'samples, fields, macro, micro, av.'], stage);
end

mvt.log(opts, 'stage %s%s finished in %.0f s', stage, dayTag(stage, day), toc(runTimer));
end

% ---------------------------------------------------------------------------
function tag = dayTag(stage, day)
if strcmp(stage, 'av') || isempty(day)
    tag = '';
else
    tag = sprintf(' 2022-11-%d', day);
end
end
