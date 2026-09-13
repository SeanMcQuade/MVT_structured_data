function report = verify_outputs(day, varargin)
% VERIFY_OUTPUTS  Check a day's outputs against its recorded manifest.
%
% Purpose
%   The full-scale counterpart to the unit suite: confirms that the pipeline
%   still produces the bytes a known-good run produced. Use it after changing a
%   stage, or to check a rebuild in a separate tree:
%
%     cd MVT_structured_data/Scripts
%     verify_outputs(16)
%     MVT_RESULTS_DIR=../results_verify verify_outputs(16)
%
% Inputs
%   day       16, 17, or 18
%   varargin  options struct and/or name/value pairs; Verbose prints each file
%
% Outputs
%   report  struct with fields checked, matched, differing, missing, and a
%           details struct array for anything that failed
%
% Notes
%   Figures are recorded as size-only and reported, never failed: pixel output
%   legitimately varies with renderer and MATLAB release.
%
% Dependencies
%   mvt.paths, mvt.options, mvt.matHash
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = mvt.options(varargin{:});
mvt.assertDay(day);
p = mvt.paths();

manifestFile = fullfile(fileparts(p.scriptsDir), 'tests', 'manifests', ...
    sprintf('manifest_2022-11-%d.json', day));
if ~isfile(manifestFile)
    error('verify_outputs:noManifest', ...
        ['No manifest for 2022-11-%d. Record one from a known-good tree:\n' ...
        '  MVT_RESULTS_DIR=../results_groundtruth generate_manifest(%d)'], day, day);
end

manifest = jsondecode(fileread(manifestFile));
entries = manifest.entries;

report = struct('checked', 0, 'matched', 0, 'differing', 0, 'missing', 0, ...
    'informational', 0, 'details', struct('file', {}, 'problem', {}));

for iEntry = 1:numel(entries)
    entry = entries(iEntry);
    file = fullfile(p.resultsDir, entry.file);
    report.checked = report.checked + 1;

    if ~isfile(file)
        report.missing = report.missing + 1;
        report.details(end+1) = struct('file', entry.file, 'problem', 'missing'); %#ok<AGROW>
        continue
    end

    switch entry.kind
        case 'json-md5'
            actual = md5OfFile(file);
        case 'mat-content'
            actual = mvt.matHash(file);
        otherwise
            % Size-only entries are informational.
            info = dir(file);
            report.informational = report.informational + 1;
            if info(1).bytes ~= entry.bytes && opts.Verbose
                fprintf('  note: %s is %d bytes, manifest recorded %d\n', ...
                    entry.file, info(1).bytes, entry.bytes);
            end
            continue
    end

    if strcmp(actual, entry.digest)
        report.matched = report.matched + 1;
        if opts.Verbose
            fprintf('  ok      %s\n', entry.file);
        end
    else
        report.differing = report.differing + 1;
        report.details(end+1) = struct('file', entry.file, 'problem', ...
            sprintf('%s mismatch (expected %s, got %s)', entry.kind, ...
            entry.digest, actual)); %#ok<AGROW>
        fprintf('  DIFFERS %s\n', entry.file);
    end
end

fprintf(['\n2022-11-%d: %d checked, %d matched, %d differing, %d missing, ' ...
    '%d informational (figures)\n'], day, report.checked, report.matched, ...
    report.differing, report.missing, report.informational);

if report.differing > 0 || report.missing > 0
    for iDetail = 1:numel(report.details)
        fprintf('  %s: %s\n', report.details(iDetail).file, report.details(iDetail).problem);
    end
    error('verify_outputs:mismatch', ...
        '%d files differ and %d are missing for 2022-11-%d', ...
        report.differing, report.missing, day);
end
end

% ---------------------------------------------------------------------------
function digest = md5OfFile(file)
if ismac
    command = sprintf('md5 -q ''%s''', file);
else
    command = sprintf('md5sum ''%s'' | cut -d" " -f1', file);
end
[status, output] = system(command);
if status ~= 0
    error('verify_outputs:md5Failed', 'Could not hash %s: %s', file, output);
end
digest = strtrim(output);
end
