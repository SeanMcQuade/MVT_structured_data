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
% tempdir rather than a literal '/tmp': macOS resolves /tmp through a symlink
% to /private/tmp, so comparing against the literal fails for a reason that has
% nothing to do with the override working.
target = fullfile(tempdir, 'mvt_results_elsewhere');
setenv('MVT_RESULTS_DIR', target);
p = mvt.paths();
verifyEqual(testCase, p.resultsDir, target);
verifyEqual(testCase, p.stateDir, fullfile(target, '.mvt'));
% Derived locations must follow the override, or a comparison run would write
% its build state into the real tree.
verifyTrue(testCase, startsWith(p.cacheDir, target));
verifyTrue(testCase, startsWith(p.stampDir, target));
end

function testRelativeOverrideResolvesAgainstWorkingDirectory(testCase)
% A relative override follows the caller's working directory, as every other
% command-line tool does, and as the Python CLI does with the same string. It
% used to resolve against the workspace root; this test asserted that older
% behaviour and outlived it.
setenv('MVT_RESULTS_DIR', 'results_verify');
p = mvt.paths();
verifyEqual(testCase, p.resultsDir, fullfile(pwd, 'results_verify'));
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
verifyEqual(testCase, outputs{1}, fullfile(p.resultsDir, 'analysis', '2022-11-16', ...
    'samples_for_distance_analysis_16.mat'));

outputs = mvt.expectedOutputs('fields', 17);
verifySubstring(testCase, outputs{1}, 'fields_motion_2022-11-17.mat');
end

function testIntermediatesAndFiguresAreSeparate(testCase)
% A results-only download takes analysis/ (the inputs) without figures/ (the
% outputs), so no stage may mix the two.
p = mvt.paths();
for stage = {'lanes', 'lc', 'samples', 'fields'}
    outputs = mvt.expectedOutputs(stage{1}, 17);
    for iOut = 1:numel(outputs)
        verifyEqual(testCase, fileparts(outputs{iOut}), ...
            fullfile(p.resultsDir, 'analysis', '2022-11-17'), ...
            sprintf('%s writes outside analysis/', stage{1}));
    end
end
for stage = {'lcplot', 'macro', 'micro'}
    outputs = mvt.expectedOutputs(stage{1}, 17);
    for iOut = 1:numel(outputs)
        verifyEqual(testCase, fileparts(outputs{iOut}), ...
            fullfile(p.resultsDir, 'figures', '2022-11-17'), ...
            sprintf('%s writes outside figures/', stage{1}));
    end
end
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
