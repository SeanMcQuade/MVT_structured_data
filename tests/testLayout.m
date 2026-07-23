function tests = testLayout
% TESTLAYOUT  Tests for path resolution, output naming, and expected outputs.
%
% These pin the conventions that the Makefile, the stages, and the tests all
% depend on agreeing about: where results live, and what a processed segment is
% called.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'Scripts'));
testCase.TestData.savedResults = getenv('MVT_RESULTS_DIR');
testCase.TestData.savedData = getenv('MVT_DATA_DIR');
end

function teardownOnce(testCase)
setenv('MVT_RESULTS_DIR', testCase.TestData.savedResults);
setenv('MVT_DATA_DIR', testCase.TestData.savedData);
end

function setup(~)
setenv('MVT_RESULTS_DIR', '');
setenv('MVT_DATA_DIR', '');
end

function teardown(~)
setenv('MVT_RESULTS_DIR', '');
setenv('MVT_DATA_DIR', '');
end

% ---------------------------------------------------------------------------
function testRepositoryIsSiblingOfDataAndResults(testCase)
p = mvt.paths();
verifyEqual(testCase, fileparts(p.repoRoot), p.dataRoot);
verifyEqual(testCase, p.resultsDir, fullfile(p.dataRoot, 'results'));
verifyEqual(testCase, p.dataDir, fullfile(p.dataRoot, 'data'));
verifyEqual(testCase, p.scriptsDir, fullfile(p.repoRoot, 'Scripts'));
end

function testResultsOverrideIsHonored(testCase)
setenv('MVT_RESULTS_DIR', '/tmp/mvt_results_elsewhere');
p = mvt.paths();
verifyEqual(testCase, p.resultsDir, '/tmp/mvt_results_elsewhere');
verifyEqual(testCase, p.stateDir, fullfile('/tmp/mvt_results_elsewhere', '.mvt'));
% Derived locations must follow the override, or a comparison run would write
% its build state into the real tree.
verifyTrue(testCase, startsWith(p.cacheDir, '/tmp/mvt_results_elsewhere'));
verifyTrue(testCase, startsWith(p.stampDir, '/tmp/mvt_results_elsewhere'));
end

function testRelativeOverrideResolvesAgainstWorkspace(testCase)
setenv('MVT_RESULTS_DIR', 'results_verify');
p = mvt.paths();
verifyEqual(testCase, p.resultsDir, fullfile(p.dataRoot, 'results_verify'));
end

function testDayFolders(testCase)
p = mvt.paths();
verifyEqual(testCase, mvt.dayDir('slim', 17), ...
    fullfile(p.resultsDir, 'slim', '2022-11-17'));
verifyEqual(testCase, mvt.dayDir('raw', 16), ...
    fullfile(p.dataDir, 'i24motion', '2022-11-16'));
verifyEqual(testCase, mvt.dayDir('gps', 18), fullfile(p.resultsDir, 'gps'));
verifyError(testCase, @() mvt.dayDir('nonsense', 16), 'mvt:dayDir:unknownKind');
end

function testSegmentNameMatchesReleasedConvention(testCase)
% 1668599999.9835877 is the first_timestamp of segment 00 on 2022-11-16, and
% the released file is named for it in Nashville local time.
name = mvt.segmentName(1668599999.9835877);
verifyEqual(testCase, name, 'I-24MOTION_2022-11-16_05-59-59.json');
end

function testSegmentNameRejectsRubbish(testCase)
verifyError(testCase, @() mvt.segmentName(NaN), 'mvt:segmentName:badTimestamp');
verifyError(testCase, @() mvt.segmentName([1 2]), 'mvt:segmentName:badTimestamp');
end

function testExpectedOutputsForDataStages(testCase)
p = mvt.paths();
outputs = mvt.expectedOutputs('samples', 16);
verifyEqual(testCase, outputs{1}, fullfile(p.resultsDir, 'figures', '2022-11-16', ...
    'samples_for_distance_analysis_16.mat'));

outputs = mvt.expectedOutputs('fields', 17);
verifySubstring(testCase, outputs{1}, 'fields_motion_2022-11-17.mat');
end

function testAvFiguresLiveInTheSharedFolder(testCase)
% They compare days, so they must not land in any one day's folder.
p = mvt.paths();
outputs = mvt.expectedOutputs('av', 16);
for iOut = 1:numel(outputs)
    verifyEqual(testCase, fileparts(outputs{iOut}), fullfile(p.resultsDir, 'figures'));
end
end

function testEnsureDirIsIdempotent(testCase)
folder = fullfile(tempname, 'nested', 'deeper');
cleanup = onCleanup(@() rmdir(fileparts(fileparts(folder)), 's')); %#ok<NASGU>
mvt.ensureDir(folder);
mvt.ensureDir(folder);   % must not raise
verifyTrue(testCase, isfolder(folder));
end

function testAtomicWriteLeavesNoTemporaryFile(testCase)
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>

target = fullfile(folder, 'out.json');
mvt.atomicWrite(target, '[{"a":1}]');

verifyEqual(testCase, fileread(target), '[{"a":1}]');
leftovers = dir(fullfile(folder, '*.tmp'));
verifyEmpty(testCase, leftovers);
end

function testAtomicWriteReplacesALinkRatherThanFollowingIt(testCase)
% Verification runs often link inputs from the real results tree. A stage that
% writes to such a path must replace the link, never write through it.
if ispc
    return   % symlinks need elevation on Windows
end
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>

protected = fullfile(folder, 'original.json');
fid = fopen(protected, 'w'); fwrite(fid, 'ORIGINAL', 'char'); fclose(fid);

link = fullfile(folder, 'link.json');
system(sprintf('ln -s ''%s'' ''%s''', protected, link));

mvt.atomicWrite(link, 'REPLACED');

verifyEqual(testCase, fileread(protected), 'ORIGINAL');
verifyEqual(testCase, fileread(link), 'REPLACED');
end
