function atomicWrite(targetPath, payload)
% MVT.ATOMICWRITE  Write a file via a temporary name, then rename into place.
%
% Purpose
%   Under `make -j` (and any interrupted run) a stage can die halfway through
%   writing a multi-hundred-megabyte JSON. The half-written file would then
%   look complete to the staleness check and poison every downstream stage.
%   Writing to a process-unique temporary name and renaming makes the target
%   appear only once it is whole.
%
% Inputs
%   targetPath  char, final path
%   payload     char row vector or uint8 vector to write verbatim
%
% Outputs
%   (none; raises mvt:atomicWrite:* on failure)
%
% Notes
%   The rename is atomic within a filesystem, which is the case here: the
%   temporary file is created in the target's own folder.
%
% Dependencies
%   mvt.ensureDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

folder = fileparts(targetPath);
mvt.ensureDir(folder);

[~, base, ext] = fileparts(targetPath);
tempPath = fullfile(folder, sprintf('.%s%s.%d.tmp', base, ext, processId()));

fid = fopen(tempPath, 'w');
if fid < 0
    error('mvt:atomicWrite:cannotOpen', ...
        'Could not open temporary file %s for writing.', tempPath);
end

try
    if ischar(payload) || isstring(payload)
        fwrite(fid, char(payload), 'char');
    else
        fwrite(fid, payload, 'uint8');
    end
    fclose(fid);
catch err
    fclose(fid);
    delete(tempPath);
    rethrow(err);
end

[ok, msg] = movefile(tempPath, targetPath, 'f');
if ~ok
    delete(tempPath);
    error('mvt:atomicWrite:cannotRename', ...
        'Could not move %s to %s: %s', tempPath, targetPath, msg);
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
