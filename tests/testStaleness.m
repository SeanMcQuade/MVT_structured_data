function tests = testStaleness
% TESTSTALENESS  Tests for mvt.isStale, the rule that decides every rebuild.
%
% Uses temporary files rather than the real results tree, so the suite runs in
% seconds on a machine with no data.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'Scripts'));
end

function setup(testCase)
testCase.TestData.folder = tempname;
mkdir(testCase.TestData.folder);
end

function teardown(testCase)
rmdir(testCase.TestData.folder, 's');
end

% ---------------------------------------------------------------------------
function testMissingOutputIsStale(testCase)
output = fullfile(testCase.TestData.folder, 'absent.json');
[stale, reason] = mvt.isStale(output, {}, {});
verifyTrue(testCase, stale);
verifySubstring(testCase, reason, 'missing output');
end

function testFreshOutputIsNotStale(testCase)
input = writeFile(testCase, 'in.json', 'a', 60);   % input is a minute old
output = writeFile(testCase, 'out.json', 'b');
[stale, reason] = mvt.isStale(output, input, {});
verifyFalse(testCase, stale);
verifySubstring(testCase, reason, 'up to date');
end

function testNewerInputIsStale(testCase)
output = writeFile(testCase, 'out.json', 'b', 60);  % output is a minute old
input = writeFile(testCase, 'in.json', 'a');
[stale, reason] = mvt.isStale(output, input, {});
verifyTrue(testCase, stale);
verifySubstring(testCase, reason, 'newer than output');
end

function testNewerSourceIsStale(testCase)
% The behavior the original pipeline lacked: editing code invalidates output.
output = writeFile(testCase, 'out.json', 'b', 60);
source = writeFile(testCase, 'stage.m', 'function y = stage(x)');
[stale, reason] = mvt.isStale(output, {}, source);
verifyTrue(testCase, stale);
verifySubstring(testCase, reason, 'stage.m');
end

function testForceOverridesFreshness(testCase)
input = writeFile(testCase, 'in.json', 'a', 60);
output = writeFile(testCase, 'out.json', 'b');
[stale, reason] = mvt.isStale(output, input, {}, mvt.options('Force', true));
verifyTrue(testCase, stale);
verifySubstring(testCase, reason, 'forced');
end

function testGlobWithNoMatchesIsStale(testCase)
pattern = fullfile(testCase.TestData.folder, 'I-24MOTION_*.json');
[stale, reason] = mvt.isStale(pattern, {}, {});
verifyTrue(testCase, stale);
verifySubstring(testCase, reason, 'no files match');
end

function testOldestOutputDecides(testCase)
% A partially rebuilt target must read as stale, not fresh.
writeFile(testCase, 'out1.json', 'x', 120);   % oldest
input = writeFile(testCase, 'in.json', 'a', 60);
writeFile(testCase, 'out2.json', 'y');       % newest
pattern = fullfile(testCase.TestData.folder, 'out*.json');
verifyTrue(testCase, mvt.isStale(pattern, input, {}));
end

function testMissingPrerequisiteIsIgnored(testCase)
% The stage itself raises the informative error about absent inputs.
output = writeFile(testCase, 'out.json', 'b');
absent = fullfile(testCase.TestData.folder, 'never_created.json');
verifyFalse(testCase, mvt.isStale(output, absent, {}));
end

function testNoDeclaredOutputsIsStale(testCase)
verifyTrue(testCase, mvt.isStale({}, {}, {}));
end

% ---------------------------------------------------------------------------
function path = writeFile(testCase, name, contents, ageSeconds)
% Write a file and, optionally, backdate it by ageSeconds.
%
% mvt.isStale compares timestamps with a one-second tolerance, and dir()
% reports them only to the second, so tests set modification times explicitly
% rather than sleeping between writes: deterministic, and instant.
if nargin < 4
    ageSeconds = 0;
end
path = fullfile(testCase.TestData.folder, name);
fid = fopen(path, 'w');
fwrite(fid, contents, 'char');
fclose(fid);

if ageSeconds ~= 0
    stamp = datestr(datetime('now') - seconds(ageSeconds), 'yyyymmddHHMM.SS'); %#ok<DATST>
    [status, message] = system(sprintf('touch -t %s ''%s''', stamp, path));
    if status ~= 0
        error('testStaleness:touchFailed', 'Could not backdate %s: %s', path, message);
    end
end
end
