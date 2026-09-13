function failed = verify(varargin)
% MVT.VERIFY  Check generated outputs against the expected checksums.
%
% Purpose
%   Answers "did this machine produce the same bytes as the reference?" from
%   inside MATLAB, so a build on a machine without Python can still be
%   validated. Reads the same manifests the Python tool uses
%   (python/expected/checksums-2022-11-DD.json), so both implementations check
%   against one set of expected values.
%
% Usage
%   mvt.verify                     % all three days
%   mvt.verify Days 18             % one day
%   mvt.verify Verbose true        % list every file, not just failures
%   failed = mvt.verify(...);      % returns the count instead of erroring
%
% Inputs
%   Name/value: 'Days' (default [16 17 18]), plus any mvt.options setting.
%
% Outputs
%   failed  number of files that did not match. Raises instead when called
%           without an output argument and anything failed.
%
% What is and is not checked
%   Only the JSON products (gps, slim) are byte-comparable and therefore
%   checked. Figures are skipped: renderers differ and PNG bytes carry encoder
%   metadata. .mat files are skipped: v7 is a gzip stream carrying a creation
%   timestamp and v7.3 is an HDF5 container, so two saves of identical data
%   differ on disk - use mvt.matHash for those. dataset_info.json is skipped
%   because it records host and user by design.
%
% Notes
%   A manifest records the data set version it was built from. If that differs
%   from mvt.dataVersion() this warns before reporting, because the mismatches
%   are then expected rather than a defect.
%
% Dependencies
%   mvt.options, mvt.paths, mvt.dayDir, mvt.dataVersion, mvt.manifest
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
p = mvt.paths();

failed = 0;
checked = 0;
matched = 0;
missing = 0;

for day = opts.Days
    manifestFile = fullfile(p.repoRoot, 'python', 'expected', ...
        sprintf('checksums-2022-11-%d.json', day));
    if ~isfile(manifestFile)
        fprintf(2, '[mvt] no checksum manifest for 2022-11-%d at %s\n', ...
            day, manifestFile);
        failed = failed + 1;
        continue
    end
    stored = jsondecode(fileread(manifestFile));

    fprintf('[mvt] verify 2022-11-%d against %s\n', day, manifestFile);
    if isfield(stored, 'data_version')
        if ~strcmp(stored.data_version, mvt.dataVersion())
            fprintf(2, ['[mvt] WARNING: manifest describes data version %s, ' ...
                'this code produces %s.\n       Differences below are ' ...
                'expected, not defects. See docs/DATA_CHANGELOG.md.\n'], ...
                stored.data_version, mvt.dataVersion());
        end
    else
        fprintf(2, ['[mvt] WARNING: manifest records no data version, so it ' ...
            'cannot be matched\n       against this code (%s).\n'], ...
            mvt.dataVersion());
    end

    names = fieldnames(stored.files);
    for iFile = 1:numel(names)
        entry = stored.files.(names{iFile});
        % The manifest carries the path explicitly: jsondecode mangles object
        % keys into struct field names ('/', '-' and '.' all become '_'),
        % which cannot be reversed.
        relative = entry.path;
        actualFile = fullfile(p.resultsDir, relative);
        checked = checked + 1;

        if ~isfile(actualFile)
            fprintf(2, '  MISSING  %s\n', relative);
            missing = missing + 1;
            failed = failed + 1;
            continue
        end
        digest = md5(actualFile);
        if strcmp(digest, entry.md5)
            matched = matched + 1;
            if opts.Verbose
                fprintf('  ok       %s\n', relative);
            end
        else
            info = dir(actualFile);
            fprintf(2, '  DIFFERS  %s\n           got %s, expected %s (%d vs %d bytes)\n', ...
                relative, digest, entry.md5, info.bytes, entry.bytes);
            failed = failed + 1;
        end
    end
end

fprintf('[mvt] %d checked: %d match, %d differ, %d missing\n', ...
    checked, matched, failed - missing, missing);
if failed == 0 && checked > 0
    fprintf('[mvt] VERIFIED\n');
end

if failed > 0 && nargout == 0
    error('mvt:verify:mismatch', '%d of %d file(s) did not match.', failed, checked);
end
end

% ---------------------------------------------------------------------------
function digest = md5(file)
% Stream the file through Java's MD5 - these run to gigabytes, so never slurp.
persistent bufferSize
if isempty(bufferSize)
    bufferSize = 1024 * 1024;
end
engine = java.security.MessageDigest.getInstance('MD5');
fid = fopen(file, 'r');
if fid < 0
    error('mvt:verify:cannotRead', 'Could not open %s', file);
end
cleanup = onCleanup(@() fclose(fid));
while true
    block = fread(fid, bufferSize, '*uint8');
    if isempty(block)
        break
    end
    engine.update(block);
end
digest = lower(reshape(dec2hex(typecast(engine.digest(), 'uint8')).', 1, []));
end

