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

banner(opts, stage, dayTag(stage, day));
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

footer(opts, stage, dayTag(stage, day), toc(runTimer));
end

% ---------------------------------------------------------------------------
function banner(opts, stage, tag)
% A visible start marker. Under `make -j` several MATLABs write to the same
% pipe, so every line is self-identifying rather than relying on position.
if ~(isstruct(opts) && isfield(opts,'Verbose') && opts.Verbose)
    return
end
label = sprintf('%s%s', stage, tag);
fprintf('\n%s\n', repmat('=', 1, 64));
fprintf('  MVT %s   started %s\n', upper(label), datestr(now, 'HH:MM:SS')); %#ok<TNOW1,DATST>
if isfield(opts,'Shard') && numel(opts.Shard) == 2 && opts.Shard(2) > 1
    fprintf('  shard %d of %d\n', opts.Shard(1), opts.Shard(2));
end
fprintf('%s\n', repmat('=', 1, 64));
end

% ---------------------------------------------------------------------------
function footer(opts, stage, tag, seconds)
if ~(isstruct(opts) && isfield(opts,'Verbose') && opts.Verbose)
    return
end
if seconds < 60
    took = sprintf('%ds', round(seconds));
elseif seconds < 3600
    took = sprintf('%dm%02ds', floor(seconds/60), round(mod(seconds,60)));
else
    took = sprintf('%dh%02dm', floor(seconds/3600), floor(mod(seconds,3600)/60));
end
fprintf('%s\n', repmat('-', 1, 64));
fprintf('  MVT %s%s   DONE in %s\n', upper(stage), upper(tag), took);
fprintf('%s\n\n', repmat('-', 1, 64));
end

% ---------------------------------------------------------------------------
function tag = dayTag(stage, day)
if strcmp(stage, 'av') || isempty(day)
    tag = '';
else
    tag = sprintf(' 2022-11-%d', day);
end
end
