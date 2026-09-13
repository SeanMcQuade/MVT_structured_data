function tests = testComparators
% TESTCOMPARATORS  Tests for mvt.matHash and mvt.matCompare.
%
% These are the tools that decide whether a rebuilt .mat output matches its
% reference, so their failure modes matter: a hash that ignores a difference,
% or a comparison that reports one where there is none, would quietly validate
% a broken pipeline.
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
function testHashIgnoresContainerDifferences(testCase)
% The reason file md5 is useless here: the same data saved twice, in two
% formats, must hash the same.
data.samples = (1:1000)' * 0.5;
data.label = 'day16';

v7 = fullfile(testCase.TestData.folder, 'v7.mat');
v73 = fullfile(testCase.TestData.folder, 'v73.mat');
save(v7, '-struct', 'data');
save(v73, '-struct', 'data', '-v7.3');

verifyEqual(testCase, mvt.matHash(v7), mvt.matHash(v73));
end

function testHashDetectsASingleChangedValue(testCase)
data.samples = ones(1000, 1);
first = fullfile(testCase.TestData.folder, 'a.mat');
save(first, '-struct', 'data');

data.samples(500) = 1 + 1e-12;
second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, '-struct', 'data');

verifyNotEqual(testCase, mvt.matHash(first), mvt.matHash(second));
end

function testHashCoversStructsAndCells(testCase)
data.field = struct('Rho', magic(4), 'Q', {{1, 'two'}});
first = fullfile(testCase.TestData.folder, 'a.mat');
save(first, '-struct', 'data');

data.field.Rho(2, 2) = 99;
second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, '-struct', 'data');

verifyNotEqual(testCase, mvt.matHash(first), mvt.matHash(second));
end

function testHashIsIndependentOfVariableOrder(testCase)
first = fullfile(testCase.TestData.folder, 'a.mat');
alpha = 1:10; beta = 11:20; %#ok<NASGU>
save(first, 'alpha', 'beta');

second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, 'beta', 'alpha');

verifyEqual(testCase, mvt.matHash(first), mvt.matHash(second));
end

function testHashRoundingOption(testCase)
data.x = 1.00000001;
first = fullfile(testCase.TestData.folder, 'a.mat');
save(first, '-struct', 'data');
data.x = 1.00000002;
second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, '-struct', 'data');

verifyNotEqual(testCase, mvt.matHash(first), mvt.matHash(second));
verifyEqual(testCase, mvt.matHash(first, 'Decimals', 6), ...
    mvt.matHash(second, 'Decimals', 6));
end

function testCompareReportsIdenticalFiles(testCase)
data.samples_dist = (1:100)';
first = fullfile(testCase.TestData.folder, 'a.mat');
second = fullfile(testCase.TestData.folder, 'b.mat');
save(first, '-struct', 'data');
save(second, '-struct', 'data');

report = mvt.matCompare(first, second, 'Verbose', false);
verifyTrue(testCase, report.identical);
verifyEmpty(testCase, report.onlyExpected);
verifyEqual(testCase, report.variables(1).maxAbsDiff, 0);
end

function testCompareLocatesTheDifference(testCase)
data.samples_dist = ones(100, 1);
first = fullfile(testCase.TestData.folder, 'a.mat');
save(first, '-struct', 'data');
data.samples_dist(42) = 1.5;
second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, '-struct', 'data');

report = mvt.matCompare(first, second, 'Verbose', false);
verifyFalse(testCase, report.identical);
verifyEqual(testCase, report.variables(1).firstDiffIndex, 42);
verifyEqual(testCase, report.variables(1).maxAbsDiff, 0.5);
verifyEqual(testCase, report.variables(1).nDiffering, 1);
end

function testCompareHonorsTolerance(testCase)
data.x = ones(10, 1);
first = fullfile(testCase.TestData.folder, 'a.mat');
save(first, '-struct', 'data');
data.x(3) = 1 + 1e-9;
second = fullfile(testCase.TestData.folder, 'b.mat');
save(second, '-struct', 'data');

verifyFalse(testCase, getfield(mvt.matCompare(first, second, 'Verbose', false), 'identical')); %#ok<GFLD>
tolerant = mvt.matCompare(first, second, 'AbsTol', 1e-6, 'RelTol', 1e-6, 'Verbose', false);
verifyTrue(testCase, tolerant.identical);
end

function testCompareTreatsNaNsAsEqual(testCase)
% Released outputs are full of NaN (null in JSON); NaN == NaN must not be
% reported as a difference.
data.x = [1; NaN; 3];
first = fullfile(testCase.TestData.folder, 'a.mat');
second = fullfile(testCase.TestData.folder, 'b.mat');
save(first, '-struct', 'data');
save(second, '-struct', 'data');

report = mvt.matCompare(first, second, 'Verbose', false);
verifyTrue(testCase, report.identical);
end

function testCompareNotesMissingVariables(testCase)
alpha = 1:10; beta = 1:10; %#ok<NASGU>
first = fullfile(testCase.TestData.folder, 'a.mat');
second = fullfile(testCase.TestData.folder, 'b.mat');
save(first, 'alpha', 'beta');
save(second, 'alpha');

report = mvt.matCompare(first, second, 'Verbose', false);
verifyFalse(testCase, report.identical);
verifyEqual(testCase, report.onlyExpected, {'beta'});
end
