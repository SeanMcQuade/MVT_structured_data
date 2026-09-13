function idx = shardIndices(n, shard)
% MVT.SHARDINDICES  Indices this worker owns when work is split across shards.
%
% Purpose
%   Deterministically split a per-segment loop across N independent MATLAB
%   processes without any coordination between them: shard k takes every Nth
%   item starting at k. Interleaving (rather than contiguous blocks) keeps the
%   work balanced when segments differ in size.
%
% Inputs
%   n      number of items in the loop (e.g. 24 segments)
%   shard  [k N] one-based shard number and shard count (default [1 1])
%
% Outputs
%   idx  row vector of indices for this shard: k, k+N, k+2N, ... <= n
%        (empty when k > n, which is legitimate for small n and large N)
%
% Dependencies
%   None (base MATLAB).
%
% Examples
%   mvt.shardIndices(24, [1 4])  ->  1 5 9 13 17 21
%   mvt.shardIndices(24, [4 4])  ->  4 8 12 16 20 24
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 2 || isempty(shard)
    shard = [1 1];
end
if ~isnumeric(shard) || numel(shard) ~= 2 || any(~isfinite(shard))
    error('mvt:shardIndices:badShard', ...
        'Shard must be a two-element numeric vector [k N].');
end

k = double(shard(1));
N = double(shard(2));
if N < 1 || k < 1 || k > N || mod(k, 1) ~= 0 || mod(N, 1) ~= 0
    error('mvt:shardIndices:badShard', ...
        'Shard [k N] requires integers with 1 <= k <= N; got [%g %g].', k, N);
end

idx = k:N:double(n);
end
