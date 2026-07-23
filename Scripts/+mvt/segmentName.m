function name = segmentName(firstTimestamp)
% MVT.SEGMENTNAME  Output filename for a processed 10-minute MOTION segment.
%
% Purpose
%   The processed JSON files are named from the first trajectory's
%   `first_timestamp` converted to Nashville local time, e.g.
%   I-24MOTION_2022-11-16_05-59-59.json. Both the stage scripts and the build
%   manifest need that name, so the conversion lives here and nowhere else -
%   if the two ever disagreed, make would rebuild forever.
%
% Inputs
%   firstTimestamp  POSIX seconds (double), from dataTemp(1).first_timestamp
%
% Outputs
%   name  char, e.g. 'I-24MOTION_2022-11-16_05-59-59.json'
%
% Notes
%   The datetime/datestr pair reproduces the original pipeline expression
%   exactly (America/Chicago, second resolution). datestr is legacy but is kept
%   deliberately: released filenames must not change.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if ~isnumeric(firstTimestamp) || ~isscalar(firstTimestamp) || ~isfinite(firstTimestamp)
    error('mvt:segmentName:badTimestamp', ...
        'firstTimestamp must be a finite numeric scalar (POSIX seconds).');
end

fileStartT = datetime(firstTimestamp, 'convertfrom', 'posixtime', ...
    'Format', 'HH:mm:ss.SSS', 'TimeZone', 'America/Chicago');
fileStartT = datestr(fileStartT, 'YYYY-mm-dd_HH-MM-SS'); %#ok<DATST>
name = ['I-24MOTION_', fileStartT, '.json'];
end
