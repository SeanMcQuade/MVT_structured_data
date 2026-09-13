function ts = rawFirstTimestamp(rawFile)
% MVT.RAWFIRSTTIMESTAMP  First trajectory's first_timestamp from a raw segment.
%
% Purpose
%   Lets the build system learn a segment's output filename without decoding
%   the whole raw file. A raw 10-minute I-24 MOTION segment is roughly 2 GB and
%   takes minutes to jsondecode; the field we need sits a few tens of kilobytes
%   into the first record.
%
% Inputs
%   rawFile  char, path to <uuid>__{wed,thu,fri}_0_NN.json
%
% Outputs
%   ts  POSIX seconds (double)
%
% Algorithm
%   1. Read a 4 MB prefix and search it for the first "first_timestamp" value.
%   2. If the first record is unusually long, retry with a 64 MB prefix.
%   3. As a last resort decode the entire file (correct, just slow) so that an
%      unexpected layout degrades in accuracy-preserving fashion.
%
% Notes
%   first_timestamp is NOT always identical to timestamp(1) - in
%   2022-11-16 segment 00 they differ in the 7th decimal - so this field must
%   be read directly rather than inferred from the timestamp array.
%
% Dependencies
%   None (base MATLAB).
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if ~isfile(rawFile)
    error('mvt:rawFirstTimestamp:missingFile', 'No such file: %s', rawFile);
end

pattern = '"first_timestamp"\s*:\s*([-+0-9.eE]+)';

for prefixBytes = [4e6, 64e6]
    chunk = readPrefix(rawFile, prefixBytes);
    token = regexp(chunk, pattern, 'tokens', 'once');
    if ~isempty(token)
        ts = str2double(token{1});
        return
    end
end

% Fallback: decode the file properly.
data = jsondecode(fileread(rawFile));
if ~isfield(data, 'first_timestamp')
    error('mvt:rawFirstTimestamp:noField', ...
        'No first_timestamp field found in %s', rawFile);
end
ts = data(1).first_timestamp;
end

% ---------------------------------------------------------------------------
function chunk = readPrefix(file, nBytes)
fid = fopen(file, 'r');
if fid < 0
    error('mvt:rawFirstTimestamp:cannotOpen', 'Could not open %s', file);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
chunk = fread(fid, nBytes, '*char')';
end
