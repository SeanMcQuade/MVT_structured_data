function requireInputs(stageName, files, opts)
% MVT.REQUIREINPUTS  Stop a stage with an explanation, not a crash.
%
% Purpose
%   Several stages cannot run on every download. The lane-change figures need
%   LC_data, which needs the lane sidecars, which are derived from the raw
%   recordings; a reader who took only the processed trajectories has none of
%   them. That is a normal, expected state, but MATLAB's own message for it is
%   "Unable to find file", which reads like a broken installation and has sent
%   people looking for a bug that is not there.
%
%   This raises an error that says what is missing, what it would have
%   produced, what is unaffected, and how to get it.
%
% Inputs
%   stageName  char, the stage asking (e.g. 'lcplot'); selects the advice below
%   files      char / cellstr, the input files the stage requires
%   opts       (optional) options struct; unused today, accepted for symmetry
%
% Outputs
%   (none; returns silently when every file is present, raises
%   mvt:missingInput otherwise)
%
% Example
%   mvt.requireInputs('lcplot', lcDataFile)
%
% Dependencies
%   mvt.paths
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 3
    opts = struct();
end %#ok<NASGU>

if ischar(files)
    files = {files};
end
missing = files(~cellfun(@isfile, files));
if isempty(missing)
    return
end

[produces, remedy] = advice(stageName);
p = mvt.paths();

lines = {};
lines{end+1} = '';
lines{end+1} = repmat('-', 1, 72);
lines{end+1} = sprintf('  Cannot build %s: an input is not in this download.', stageName);
lines{end+1} = repmat('-', 1, 72);
lines{end+1} = '';
lines{end+1} = '  Missing:';
for iFile = 1:numel(missing)
    lines{end+1} = sprintf('    %s', shorten(missing{iFile}, p)); %#ok<AGROW>
end
lines{end+1} = '';
lines{end+1} = sprintf('  Without it this stage cannot produce %s.', produces);
lines{end+1} = '';
lines{end+1} = '  Nothing is broken, and nothing else is affected: every other';
lines{end+1} = '  stage is independent of this one, and the figures that do not';
lines{end+1} = '  depend on the missing file are unaffected. Use `make -k` to let';
lines{end+1} = '  the rest of a run finish.';
lines{end+1} = '';
lines{end+1} = '  To get it:';
for iLine = 1:numel(remedy)
    lines{end+1} = sprintf('    %s', remedy{iLine}); %#ok<AGROW>
end
lines{end+1} = '';
lines{end+1} = '  See the download routes in README.md.';
lines{end+1} = repmat('-', 1, 72);

error('mvt:missingInput', '%s', strjoin(lines, newline));
end

% ---------------------------------------------------------------------------
function [produces, remedy] = advice(stageName)
% What this stage makes, and the cheapest ways to supply what it needs.
switch lower(char(stageName))
    case 'lcplot'
        produces = 'the four lane-change figures for this day';
        remedy = { ...
            'download LC_data_DD.mat (40 MB for all three days), or', ...
            'download the lane sidecars (*_orig_dist_lane.mat, 3.6 MB) and', ...
            '  run `make lc-DD` to rebuild LC_data from them, or', ...
            'download the raw data and run `make lanes-DD lc-DD`.'};
    case 'lc'
        produces = 'LC_data_DD.mat, and so the lane-change figures';
        remedy = { ...
            'download the lane sidecars (*_orig_dist_lane.mat, 3.6 MB for', ...
            '  all three days) into results/analysis/2022-11-DD/, or', ...
            'download LC_data_DD.mat directly (40 MB) and skip this stage, or', ...
            'download the raw data and run `make lanes-DD`.'};
    case 'lanes'
        produces = 'the lane origin/destination sidecars';
        remedy = { ...
            'this stage reads the raw I-24 MOTION recordings, which are the', ...
            '  largest download. If you only want the lane-change figures,', ...
            '  download the sidecars (3.6 MB) or LC_data_DD.mat (40 MB)', ...
            '  instead of rebuilding them.'};
    case 'relspeedplot'
        produces = 'the two relative-speed figures for this day';
        remedy = { ...
            'download relspeed_data_DD.mat (690 MB for all three days), or', ...
            'download the processed trajectories (results/slim/) and run', ...
            '  `make relspeed-DD` to rebuild it.'};
    otherwise
        produces = 'its outputs';
        remedy = {'see the download routes in README.md.'};
end
end

% ---------------------------------------------------------------------------
function short = shorten(file, p)
% Paths relative to the workspace read better in a message than absolute ones.
short = file;
prefix = [p.dataRoot, filesep];
if startsWith(short, prefix)
    short = short(numel(prefix)+1:end);
end
end
