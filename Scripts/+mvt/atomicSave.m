function atomicSave(targetPath, vars, varargin)
% MVT.ATOMICSAVE  Write a .mat file through a temporary name and rename.
%
% Purpose
%   The .mat counterpart of mvt.atomicWrite. A stage killed part way through
%   `save` leaves a truncated file that still satisfies `isfile`, so the
%   staleness check calls it fresh and every downstream stage reads garbage.
%   Writing to a process-unique temporary name in the target's own folder and
%   renaming means the target appears only once it is whole.
%
% Inputs
%   targetPath  char, final path (normally ending in .mat)
%   vars        scalar struct; each field is saved as a variable of that name,
%               exactly as `save('-struct', ...)` does
%   ...         (optional) further arguments passed to `save`, e.g. '-v7.3'
%
% Outputs
%   (none; raises mvt:atomicSave:* on failure)
%
% Notes
%   The temporary name starts with a dot and ends in .tmp so that neither
%   dir('*.mat') nor the dataset_info writers can mistake a half-written file
%   for a product. The rename is atomic within a filesystem, which is the case
%   here because the temporary file is created alongside its target.
%
% Dependencies
%   mvt.ensureDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if ~isstruct(vars) || ~isscalar(vars)
    error('mvt:atomicSave:badVars', ...
        'vars must be a scalar struct whose fields become saved variables.');
end

folder = fileparts(targetPath);
if isempty(folder)
    error('mvt:atomicSave:relativePath', ...
        'targetPath must include a folder: %s', targetPath);
end
mvt.ensureDir(folder);

[~, base, ext] = fileparts(targetPath);
tempPath = fullfile(folder, sprintf('.%s%s.%d.tmp', base, ext, processId()));

try
    save(tempPath, '-struct', 'vars', varargin{:});
catch err
    deleteIfPresent(tempPath);
    rethrow(err);
end

[ok, msg] = movefile(tempPath, targetPath, 'f');
if ~ok
    deleteIfPresent(tempPath);
    error('mvt:atomicSave:cannotRename', ...
        'Could not move %s to %s: %s', tempPath, targetPath, msg);
end
end

% ---------------------------------------------------------------------------
function deleteIfPresent(file)
if isfile(file)
    delete(file);
end
end

% ---------------------------------------------------------------------------
function pid = processId()
try
    pid = feature('getpid');
catch
    % feature() is undocumented; fall back to a random-ish unique suffix.
    [~, name] = fileparts(tempname);
    pid = abs(sum(double(name)));
end
end
