function tests = testOptions
% TESTOPTIONS  Tests for mvt.options and mvt.shardIndices.
%
% These cover the switches that decide whether work happens at all, so a
% regression here silently changes what the pipeline rebuilds.
%
% Run with:  runtests('tests')
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'Scripts'));
testCase.TestData.savedEnv = saveEnvironment();
end

function teardownOnce(testCase)
restoreEnvironment(testCase.TestData.savedEnv);
end

function setup(~)
clearEnvironment();
end

function teardown(~)
clearEnvironment();
end

% ---------------------------------------------------------------------------
function testDefaults(testCase)
opts = mvt.options();
verifyFalse(testCase, opts.Force);
verifyFalse(testCase, opts.Clean);
verifyFalse(testCase, opts.DryRun);
verifyTrue(testCase, opts.Verbose);
verifyEqual(testCase, opts.Shard, [1 1]);
verifyEqual(testCase, opts.Days, [16 17 18]);
verifyEqual(testCase, opts.SettleSeconds, 5);
end

function testNameValuePairs(testCase)
opts = mvt.options('Force', true, 'Days', 17, 'SettleSeconds', 0);
verifyTrue(testCase, opts.Force);
verifyEqual(testCase, opts.Days, 17);
verifyEqual(testCase, opts.SettleSeconds, 0);
end

function testStructPassesThrough(testCase)
% Stages receive an options struct from mvt.build and re-normalize it; that
% round trip must be lossless.
first = mvt.options('Force', true, 'Shard', [2 4]);
second = mvt.options(first);
verifyEqual(testCase, second.Force, true);
verifyEqual(testCase, second.Shard, [2 4]);
end

function testUnknownOptionIsRejected(testCase)
% A typo must not silently disable a rebuild.
verifyError(testCase, @() mvt.options('Forced', true), 'mvt:options:unknownOption');
end

function testEnvironmentOverrides(testCase)
setenv('MVT_FORCE', '1');
setenv('MVT_SHARD', '3/4');
setenv('MVT_DAYS', '16 18');
opts = mvt.options();
verifyTrue(testCase, opts.Force);
verifyEqual(testCase, opts.Shard, [3 4]);
verifyEqual(testCase, opts.Days, [16 18]);
end

function testExplicitArgumentsBeatEnvironment(testCase)
setenv('MVT_FORCE', '1');
opts = mvt.options('Force', false);
verifyFalse(testCase, opts.Force);
end

function testBadEnvironmentValueIsRejected(testCase)
setenv('MVT_FORCE', 'perhaps');
verifyError(testCase, @() mvt.options(), 'mvt:options:badEnvFlag');
end

function testShardIndicesInterleave(testCase)
verifyEqual(testCase, mvt.shardIndices(24, [1 4]), 1:4:24);
verifyEqual(testCase, mvt.shardIndices(24, [4 4]), 4:4:24);
verifyEqual(testCase, mvt.shardIndices(24, [1 1]), 1:24);
end

function testShardIndicesCoverEveryItemExactlyOnce(testCase)
% The correctness property that makes sharded runs equal serial runs.
covered = [];
for k = 1:5
    covered = [covered, mvt.shardIndices(24, [k 5])]; %#ok<AGROW>
end
verifyEqual(testCase, sort(covered), 1:24);
end

function testShardIndicesRejectsBadSpec(testCase)
verifyError(testCase, @() mvt.shardIndices(10, [5 4]), 'mvt:shardIndices:badShard');
verifyError(testCase, @() mvt.shardIndices(10, [0 4]), 'mvt:shardIndices:badShard');
end

function testAssertDay(testCase)
verifyError(testCase, @() mvt.assertDay(15), 'mvt:assertDay:badDay');
verifyError(testCase, @() mvt.assertDay('16'), 'mvt:assertDay:badDay');
mvt.assertDay(16);   % must not raise
end

% ---------------------------------------------------------------------------
function saved = saveEnvironment()
names = {'MVT_FORCE', 'MVT_CLEAN', 'MVT_DRYRUN', 'MVT_VERBOSE', 'MVT_SHARD', ...
    'MVT_DAYS', 'MVT_SETTLE_SECONDS'};
saved = struct('names', {names}, 'values', {cellfun(@getenv, names, 'UniformOutput', false)});
end

function restoreEnvironment(saved)
for iName = 1:numel(saved.names)
    setenv(saved.names{iName}, saved.values{iName});
end
end

function clearEnvironment()
for name = {'MVT_FORCE', 'MVT_CLEAN', 'MVT_DRYRUN', 'MVT_VERBOSE', 'MVT_SHARD', ...
        'MVT_DAYS', 'MVT_SETTLE_SECONDS'}
    setenv(name{1}, '');
end
end
