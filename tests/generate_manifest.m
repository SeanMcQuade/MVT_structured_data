function manifestFile = generate_manifest(day, varargin)
% GENERATE_MANIFEST  Record checksums of a day's outputs for later verification.
%
% Purpose
%   Captures what a known-good run produced, so a later run can be checked
%   against it without keeping a second copy of the data. Each output kind gets
%   the comparison that suits it:
%     JSON    md5 of the bytes - these are deterministic
%     .mat    mvt.matHash content hash - file bytes are not comparable, since
%             v7 files carry a gzip timestamp and v7.3 files are HDF5
%     figures size only, recorded for information; pixels legitimately vary
%             with renderer and MATLAB release
%
% Inputs
%   day       16, 17, or 18
%   varargin  options struct and/or name/value pairs (see mvt.options)
%
% Outputs
%   manifestFile  path to tests/manifests/manifest_2022-11-DD.json
%
% Usage
%   cd MVT_structured_data/Scripts
%   generate_manifest(16)                              % from the current results
%   MVT_RESULTS_DIR=../results_groundtruth generate_manifest(16)
%
% Dependencies
%   mvt.paths, mvt.options, mvt.expectedOutputs, mvt.matHash, mvt.ensureDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
mvt.assertDay(day);
p = mvt.paths();

stages = {'gps', 'slim', 'full', 'samples', 'fields', 'macro', 'micro'};
entries = struct('stage', {}, 'file', {}, 'kind', {}, 'bytes', {}, 'digest', {});

for iStage = 1:numel(stages)
    outputs = mvt.expectedOutputs(stages{iStage}, day, opts);
    files = expandFiles(outputs);
    for iFile = 1:numel(files)
        entries(end+1) = describe(stages{iStage}, files{iFile}, p); %#ok<AGROW>
    end
end

manifestDir = fullfile(fileparts(p.scriptsDir), 'tests', 'manifests');
mvt.ensureDir(manifestDir);
manifestFile = fullfile(manifestDir, sprintf('manifest_2022-11-%d.json', day));

payload = struct('day', day, 'release', version('-release'), ...
    'resultsDir', p.resultsDir, 'entries', entries);
fid = fopen(manifestFile, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(payload), 'char');

fprintf('Recorded %d files for 2022-11-%d in %s\n', numel(entries), day, manifestFile);
end

% ---------------------------------------------------------------------------
function entry = describe(stage, file, p)
info = dir(file);
[~, ~, ext] = fileparts(file);

switch lower(ext)
    case '.json'
        kind = 'json-md5';
        digest = md5OfFile(file);
    case '.mat'
        kind = 'mat-content';
        digest = mvt.matHash(file);
    otherwise
        kind = 'size-only';
        digest = '';
end

entry = struct('stage', stage, 'file', relativePath(file, p), 'kind', kind, ...
    'bytes', info(1).bytes, 'digest', digest);
fprintf('  %-8s %-10s %s\n', stage, kind, entry.file);
end

% ---------------------------------------------------------------------------
function files = expandFiles(specs)
files = {};
for iSpec = 1:numel(specs)
    spec = specs{iSpec};
    folder = fileparts(spec);
    listing = dir(spec);
    listing = listing(~[listing.isdir]);
    for iFile = 1:numel(listing)
        files{end+1} = fullfile(folder, listing(iFile).name); %#ok<AGROW>
    end
end
files = sort(files);
end

% ---------------------------------------------------------------------------
function short = relativePath(file, p)
prefix = [p.resultsDir, filesep];
short = file;
if startsWith(short, prefix)
    short = short(numel(prefix)+1:end);
end
end

% ---------------------------------------------------------------------------
function digest = md5OfFile(file)
% Shell out: MATLAB has no public file-hashing function, and md5/md5sum reads
% the file in one streaming pass rather than loading gigabytes into memory.
if ismac
    command = sprintf('md5 -q ''%s''', file);
else
    command = sprintf('md5sum ''%s'' | cut -d" " -f1', file);
end
[status, output] = system(command);
if status ~= 0
    error('generate_manifest:md5Failed', 'Could not hash %s: %s', file, output);
end
digest = strtrim(output);
end
