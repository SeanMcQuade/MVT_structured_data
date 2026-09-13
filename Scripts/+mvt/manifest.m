function segments = manifest(day, opts)
% MVT.MANIFEST  Map each raw MOTION segment to the processed file it produces.
%
% Purpose
%   The processed filenames are content-derived: a stage only learns that
%   `..._wed_0_07.json` becomes `I-24MOTION_2022-11-16_07-09-59.json` after
%   decoding the raw file. Make needs that mapping *before* running anything,
%   and shards need it to skip work that is already done. This builds the
%   mapping cheaply (a prefix scan per raw file) and caches it.
%
% Inputs
%   day   16, 17, or 18
%   opts  (optional) options struct; Force rebuilds the cache, Verbose logs
%
% Outputs
%   segments  1xN struct array, sorted by raw filename, with fields:
%               seq             segment index parsed from the raw name (0..23)
%               rawName         raw filename
%               rawPath         absolute path to the raw file
%               firstTimestamp  POSIX seconds of the first trajectory
%               outputName      processed filename (same for slim and full)
%
% Algorithm
%   1. List <data>/i24motion/2022-11-DD/*_{wed,thu,fri}_0_*.json, matching the
%      glob the stage scripts use.
%   2. Reuse results/.mvt/manifests/segments_2022-11-DD.json when the raw
%      listing (names, sizes, timestamps) is unchanged.
%   3. Otherwise read each file's first_timestamp (mvt.rawFirstTimestamp) and
%      convert it with mvt.segmentName - the same call the stages make.
%   4. Cache the mapping; caching failures are non-fatal.
%
% Dependencies
%   mvt.paths, mvt.dayDir, mvt.options, mvt.rawFirstTimestamp,
%   mvt.segmentName, mvt.ensureDir, mvt.log
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 2 || isempty(opts)
    opts = mvt.options();
end
mvt.assertDay(day);

p = mvt.paths();
rawDir = mvt.dayDir('raw', day);
if ~isfolder(rawDir)
    error('mvt:manifest:missingRawFolder', ...
        'Raw MOTION folder does not exist: %s', rawDir);
end

dayAbbrvs = {'wed', 'thu', 'fri'};
dayAbbrv = dayAbbrvs{day - 15};
listing = dir(fullfile(rawDir, ['*_', dayAbbrv, '_0_*.json']));
listing = listing(~[listing.isdir]);
% Skip macOS AppleDouble/metadata siblings such as ._<uuid>__wed_0_00.json,
% which match the glob but are not data (the original stage scripts filter
% these too).
listing = listing(~startsWith({listing.name}, '.'));
if isempty(listing)
    error('mvt:manifest:noRawFiles', ...
        'No raw MOTION files matching *_%s_0_*.json in %s', dayAbbrv, rawDir);
end
[~, order] = sort({listing.name});
listing = listing(order);

signature = listingSignature(listing);
cacheFile = fullfile(p.manifestDir, sprintf('segments_2022-11-%d.json', day));

if ~opts.Force
    cached = readCache(cacheFile, signature);
    if ~isempty(cached)
        segments = cached;
        return
    end
end

mvt.log(opts, 'building segment manifest for 2022-11-%d (%d raw files)', ...
    day, numel(listing));

segments = struct('seq', {}, 'rawName', {}, 'rawPath', {}, ...
    'firstTimestamp', {}, 'outputName', {});
for iFile = 1:numel(listing)
    rawPath = fullfile(rawDir, listing(iFile).name);
    ts = mvt.rawFirstTimestamp(rawPath);
    seqToken = regexp(listing(iFile).name, '_0_(\d+)\.json$', 'tokens', 'once');
    if isempty(seqToken)
        seq = iFile - 1;
    else
        seq = str2double(seqToken{1});
    end
    segments(iFile) = struct( ...
        'seq', seq, ...
        'rawName', listing(iFile).name, ...
        'rawPath', rawPath, ...
        'firstTimestamp', ts, ...
        'outputName', mvt.segmentName(ts)); %#ok<AGROW>
end

writeCache(cacheFile, day, signature, segments, opts);
end

% ---------------------------------------------------------------------------
function sig = listingSignature(listing)
% Cheap change detector for the raw folder: names, sizes, and mtimes.
parts = cell(1, numel(listing));
for iFile = 1:numel(listing)
    parts{iFile} = sprintf('%s|%d|%.6f', listing(iFile).name, ...
        listing(iFile).bytes, listing(iFile).datenum);
end
sig = strjoin(parts, ';');
end

% ---------------------------------------------------------------------------
function segments = readCache(cacheFile, signature)
segments = [];
if ~isfile(cacheFile)
    return
end
try
    cached = jsondecode(fileread(cacheFile));
catch
    return
end
if ~isfield(cached, 'signature') || ~strcmp(cached.signature, signature)
    return
end
if ~isfield(cached, 'segments') || isempty(cached.segments)
    return
end
% jsondecode returns a struct array for homogeneous records.
segments = reshape(cached.segments, 1, []);
end

% ---------------------------------------------------------------------------
function writeCache(cacheFile, day, signature, segments, opts)
try
    mvt.ensureDir(fileparts(cacheFile));
    payload = struct( ...
        'day', day, ...
        'signature', signature, ...
        'release', version('-release'), ...
        'segments', segments);
    fid = fopen(cacheFile, 'w');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid, jsonencode(payload), 'char');
catch err
    mvt.log(opts, 'could not cache the segment manifest (%s)', err.message);
end
end
