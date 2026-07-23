function log(opts, fmt, varargin)
% MVT.LOG  Print a pipeline progress line, honoring opts.Verbose.
%
% Purpose
%   Uniform, greppable progress output across stages and across the sharded
%   MATLAB processes that the Makefile launches. Each line is tagged with the
%   shard when the run is sharded, so interleaved output stays readable.
%
% Inputs
%   opts  options struct from mvt.options (Verbose, Shard, DryRun)
%   fmt   printf-style format string
%   ...   format arguments
%
% Outputs
%   (none; writes to stdout)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 1 || ~isstruct(opts) || ~isfield(opts, 'Verbose') || ~opts.Verbose
    return
end

tag = '[mvt]';
if isfield(opts, 'Shard') && numel(opts.Shard) == 2 && opts.Shard(2) > 1
    tag = sprintf('[mvt %d/%d]', opts.Shard(1), opts.Shard(2));
end
if isfield(opts, 'DryRun') && opts.DryRun
    tag = [tag ' (dry-run)'];
end

fprintf('%s %s\n', tag, sprintf(fmt, varargin{:}));
end
