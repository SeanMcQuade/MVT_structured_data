function infoFile = writeDatasetInfo(stage, day, opts)
% MVT.WRITEDATASETINFO  Write the dataset_info.json sidecar for a product.
%
% Purpose
%   Records, next to the data itself, which version of the data set a folder
%   holds and how it was produced. Two copies of `slim` at 2.1 and 2.1.1 differ
%   only in the last decimal of a few fuel totals and are otherwise
%   indistinguishable by inspection, so without this a recipient cannot tell
%   them apart.
%
% Inputs
%   stage  'gps' | 'slim' | 'full' | 'samples' | 'fields' | 'macro' | 'micro'
%   day    16, 17, or 18 (ignored for stages whose product is not per-day)
%   opts   (optional) options struct; Verbose is honored, DryRun suppresses
%
% Outputs
%   infoFile  absolute path written, or '' when nothing was written
%
% Notes
%   `generated_utc` is taken from the newest output in the folder, not from the
%   clock, so re-running this on an unchanged product does not invent a new
%   timestamp.
%
%   The file mixes reproducible facts (version, product, copyright) with
%   run provenance (host, user, MATLAB release), so it is deliberately NOT
%   byte-reproducible across machines and is excluded from the checksum
%   manifests - see mvtpy.verify, which skips it by name and says so.
%
% Dependencies
%   mvt.paths, mvt.dayDir, mvt.dataVersion, mvt.options, mvt.ensureDir
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 3 || isempty(opts)
    opts = mvt.options();
end
infoFile = '';
if isfield(opts, 'DryRun') && opts.DryRun
    return
end

[productDir, product] = productFolder(stage, day);
if isempty(productDir) || ~isfolder(productDir)
    return
end

listing = dir(fullfile(productDir, '*'));
listing = listing(~[listing.isdir] & ~startsWith({listing.name}, '.'));
listing = listing(~strcmp({listing.name}, 'dataset_info.json'));
if isempty(listing)
    return          % nothing produced here yet; do not claim otherwise
end
newest = max([listing.datenum]);

info = struct( ...
    'dataset', 'CIRCLES MegaVanderTest (MVT) derived data', ...
    'product', product, ...
    'data_version', mvt.dataVersion(), ...
    'version_scheme', ['MAJOR.MINOR.PATCH; PATCH keeps every field, meaning ' ...
                       'and format identical and changes only the last decimal ' ...
                       'of a small number of values. See docs/DATA_CHANGELOG.md.'], ...
    'upstream_motion_data', ['I-24 MOTION, versioned independently by the ' ...
                             'observatory and unmodified by this pipeline'], ...
    'copyright', '(C) 2026 CIRCLES Consortium', ...
    'license', 'BSD-3-Clause', ...
    'documentation', 'docs/DATA_CHANGELOG.md, docs/REPRODUCIBLE_QUADRATURE.md');

if ~isempty(day)
    info.day = sprintf('2022-11-%d', day);
end
info.files = numel(listing);
info.generated_utc = char(datetime(newest, 'ConvertFrom', 'datenum', ...
    'TimeZone', 'local', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''', ...
    'TimeZone', 'UTC'));

info.generated_by = runProvenance();

target = fullfile(productDir, 'dataset_info.json');
try
    mvt.ensureDir(productDir);
    fid = fopen(target, 'w');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, jsonencode(info, 'PrettyPrint', true), 'char');
    infoFile = target;
    mvt.log(opts, 'wrote %s (data version %s)', target, info.data_version);
catch err
    mvt.log(opts, 'could not write %s (%s)', target, err.identifier);
end
end

% ---------------------------------------------------------------------------
function provenance = runProvenance()
% Who and what produced this. Best effort: a missing field must never fail a run.
provenance = struct( ...
    'implementation', 'MATLAB', ...
    'matlab_release', version('-release'), ...
    'matlab_version', version(), ...
    'platform', computer(), ...
    'host', '', 'user', '', 'code_commit', '');
try
    if ispc
        provenance.host = strtrim(getenv('COMPUTERNAME'));
        provenance.user = strtrim(getenv('USERNAME'));
    else
        [~, hostOut] = system('hostname');
        provenance.host = strtrim(hostOut);
        provenance.user = strtrim(getenv('USER'));
    end
catch
end
try
    p = mvt.paths();
    [status, sha] = system(sprintf('git -C "%s" rev-parse --short HEAD', p.repoRoot));
    if status == 0
        provenance.code_commit = strtrim(sha);
    end
catch
end
end

% ---------------------------------------------------------------------------
function [folder, product] = productFolder(stage, day)
% Where a stage's released output lives, and what to call that product.
folder = '';
product = stage;
switch lower(char(stage))
    case 'gps'
        p = mvt.paths();
        folder = fullfile(p.resultsDir, 'gps');
        product = 'gps (control-vehicle 10 Hz GPS)';
    case {'slim', 'full'}
        folder = mvt.dayDir(lower(char(stage)), day);
        product = sprintf('%s (I-24 MOTION trajectories)', lower(char(stage)));
    case {'samples', 'fields', 'macro', 'micro', 'av'}
        folder = mvt.dayDir('figures', day);
        product = 'figures and analysis products';
end
end
