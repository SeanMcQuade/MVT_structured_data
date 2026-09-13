function opts = options(varargin)
% MVT.OPTIONS  Normalize pipeline run options (staleness, sharding, verbosity).
%
% Purpose
%   Give every stage function one consistent, validated options struct, and
%   let the Makefile drive the same switches through environment variables.
%
% Inputs
%   Accepts, in any combination:
%     - an options struct (typically one returned by a previous call), and/or
%     - name/value pairs.
%   Field names are matched case-insensitively; unknown names raise an error so
%   that typos such as 'Forced' cannot silently disable a rebuild.
%
% Outputs
%   opts  struct with fields:
%     Force      logical  rebuild even when outputs look up to date  (false)
%     Clean      logical  delete this stage's outputs before running (false)
%     DryRun     logical  report what would be done, write nothing   (false)
%     Verbose    logical  print progress                            (true)
%     Shard      [k N]    process shard k of N (see mvt.shardIndices) ([1 1])
%     UseParfor  logical  use parfor when the toolbox is licensed    (false)
%     Days       vector   days a driver should process           ([16 17 18])
%     SettleSeconds  numeric  pause after clearing large variables, giving
%                    the OS time to reclaim memory before the next big
%                    allocation; 0 disables                              (5)
%
% Environment overrides (applied first, so explicit arguments win)
%   MVT_FORCE, MVT_CLEAN, MVT_DRYRUN  '1'/'true'/'yes' enable
%   MVT_VERBOSE                       '0'/'false'/'no' disable
%   MVT_SHARD                         'k/N' or 'k,N'
%   MVT_DAYS                          e.g. '16 17 18' or '17'
%   MVT_SETTLE_SECONDS                e.g. '0' to disable the memory pause
%
% Algorithm
%   1. Start from the documented defaults.
%   2. Apply environment overrides (Makefile / shell driven).
%   3. Apply struct arguments, then name/value pairs, in the order given.
%   4. Validate types and the shard specification.
%
% Dependencies
%   mvt.shardIndices (validation only)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

opts = struct( ...
    'Force', false, ...
    'Clean', false, ...
    'DryRun', false, ...
    'Verbose', true, ...
    'Shard', [1 1], ...
    'UseParfor', false, ...
    'Days', [16 17 18], ...
    'SettleSeconds', 5);

known = fieldnames(opts);

% ---- environment overrides -------------------------------------------------
opts.Force   = envFlag('MVT_FORCE',   opts.Force);
opts.Clean   = envFlag('MVT_CLEAN',   opts.Clean);
opts.DryRun  = envFlag('MVT_DRYRUN',  opts.DryRun);
opts.Verbose = envFlag('MVT_VERBOSE', opts.Verbose);

shardEnv = strtrim(getenv('MVT_SHARD'));
if ~isempty(shardEnv)
    parts = str2double(regexp(shardEnv, '[/,: ]+', 'split'));
    if numel(parts) ~= 2 || any(isnan(parts))
        error('mvt:options:badEnvShard', ...
            'MVT_SHARD must look like ''k/N''; got ''%s''.', shardEnv);
    end
    opts.Shard = parts(:)';
end

settleEnv = strtrim(getenv('MVT_SETTLE_SECONDS'));
if ~isempty(settleEnv)
    settle = str2double(settleEnv);
    if isnan(settle)
        error('mvt:options:badEnvSettle', ...
            'MVT_SETTLE_SECONDS must be numeric; got ''%s''.', settleEnv);
    end
    opts.SettleSeconds = settle;
end

daysEnv = strtrim(getenv('MVT_DAYS'));
if ~isempty(daysEnv)
    parsed = str2double(regexp(daysEnv, '[,; ]+', 'split'));
    if isempty(parsed) || any(isnan(parsed))
        error('mvt:options:badEnvDays', ...
            'MVT_DAYS must be a list of days; got ''%s''.', daysEnv);
    end
    opts.Days = parsed(:)';
end

% ---- explicit arguments ----------------------------------------------------
args = varargin;
while ~isempty(args)
    if isstruct(args{1})
        given = args{1};
        args(1) = [];
        names = fieldnames(given);
        for iName = 1:numel(names)
            opts = assign(opts, known, names{iName}, given.(names{iName}));
        end
    elseif (ischar(args{1}) || isstring(args{1})) && numel(args) >= 2
        opts = assign(opts, known, char(args{1}), args{2});
        args(1:2) = [];
    elseif isempty(args{1})
        args(1) = [];   % tolerate [] placeholders from callers
    else
        error('mvt:options:badArgs', ...
            ['Options must be given as a struct or as name/value pairs; ', ...
            'received a %s with no value.'], class(args{1}));
    end
end

% ---- validation ------------------------------------------------------------
for iFlag = {'Force', 'Clean', 'DryRun', 'Verbose', 'UseParfor'}
    name = iFlag{1};
    value = opts.(name);
    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
        error('mvt:options:badFlag', '%s must be a logical scalar.', name);
    end
    opts.(name) = logical(value);
end

if ~isnumeric(opts.SettleSeconds) || ~isscalar(opts.SettleSeconds) || ...
        ~isfinite(opts.SettleSeconds) || opts.SettleSeconds < 0
    error('mvt:options:badSettle', ...
        'SettleSeconds must be a finite, non-negative numeric scalar.');
end
opts.SettleSeconds = double(opts.SettleSeconds);

if ~isnumeric(opts.Shard) || numel(opts.Shard) ~= 2
    error('mvt:options:badShard', 'Shard must be a two-element vector [k N].');
end
opts.Shard = double(opts.Shard(:)');
mvt.shardIndices(0, opts.Shard);   % validates k and N

opts.Days = double(opts.Days(:)');
end

% ---------------------------------------------------------------------------
function opts = assign(opts, known, name, value)
match = known(strcmpi(known, name));
if isempty(match)
    error('mvt:options:unknownOption', ...
        'Unknown option ''%s''. Valid options: %s.', name, strjoin(known', ', '));
end
opts.(match{1}) = value;
end

% ---------------------------------------------------------------------------
function value = envFlag(name, value)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    return
end
switch raw
    case {'1', 'true', 'yes', 'on'}
        value = true;
    case {'0', 'false', 'no', 'off'}
        value = false;
    otherwise
        error('mvt:options:badEnvFlag', ...
            '%s must be 1/0, true/false, yes/no, or on/off; got ''%s''.', ...
            name, raw);
end
end
