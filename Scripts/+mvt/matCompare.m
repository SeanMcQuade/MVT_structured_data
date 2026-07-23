function report = matCompare(expectedFile, actualFile, varargin)
% MVT.MATCOMPARE  Compare two .mat outputs variable by variable.
%
% Purpose
%   When mvt.matHash reports a difference, this says *what* differs: which
%   variables, by how much, and at which index. That is what turns "the samples
%   file changed" into a debuggable statement.
%
% Inputs
%   expectedFile  reference .mat (e.g. from results_groundtruth)
%   actualFile    .mat produced by the current code
%   varargin      name/value pairs:
%                   'AbsTol'   absolute tolerance for numeric data (default 0)
%                   'RelTol'   relative tolerance (default 0)
%                   'Verbose'  print a summary table (default true)
%
% Outputs
%   report  struct with fields:
%             identical      logical, true when nothing differs beyond tolerance
%             onlyExpected   cellstr of variables missing from the actual file
%             onlyActual     cellstr of variables the reference does not have
%             variables      struct array, one entry per compared variable:
%               name, class, sizeMatch, equal, maxAbsDiff, maxRelDiff,
%               firstDiffIndex, nDiffering, note
%
% Usage
%   report = mvt.matCompare( ...
%       '../results_groundtruth/figures/2022-11-16/samples_for_distance_analysis_16.mat', ...
%       '../results_verify/figures/2022-11-16/samples_for_distance_analysis_16.mat');
%
% Dependencies
%   mvt.log (for printing)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

options = parseOptions(varargin{:});

expected = load(expectedFile);
actual = load(actualFile);

expectedNames = sort(fieldnames(expected));
actualNames = sort(fieldnames(actual));
shared = intersect(expectedNames, actualNames);

report = struct();
report.onlyExpected = setdiff(expectedNames, actualNames);
report.onlyActual = setdiff(actualNames, expectedNames);
report.variables = struct('name', {}, 'class', {}, 'sizeMatch', {}, 'equal', {}, ...
    'maxAbsDiff', {}, 'maxRelDiff', {}, 'firstDiffIndex', {}, 'nDiffering', {}, ...
    'note', {});

for iName = 1:numel(shared)
    name = shared{iName};
    report.variables(end+1) = compareOne(name, expected.(name), actual.(name), options); %#ok<AGROW>
end

report.identical = isempty(report.onlyExpected) && isempty(report.onlyActual) && ...
    all([report.variables.equal]);

if options.Verbose
    printReport(report, expectedFile, actualFile);
end
end

% ---------------------------------------------------------------------------
function entry = compareOne(name, expectedValue, actualValue, options)
entry = struct('name', name, 'class', class(expectedValue), 'sizeMatch', true, ...
    'equal', false, 'maxAbsDiff', NaN, 'maxRelDiff', NaN, ...
    'firstDiffIndex', NaN, 'nDiffering', NaN, 'note', '');

if ~strcmp(class(expectedValue), class(actualValue))
    entry.note = sprintf('class differs: %s vs %s', class(expectedValue), class(actualValue));
    return
end
if ~isequal(size(expectedValue), size(actualValue))
    entry.sizeMatch = false;
    entry.note = sprintf('size differs: %s vs %s', ...
        mat2str(size(expectedValue)), mat2str(size(actualValue)));
    return
end

if isnumeric(expectedValue) || islogical(expectedValue)
    a = double(expectedValue(:));
    b = double(actualValue(:));
    bothNaN = isnan(a) & isnan(b);
    differing = ~(a == b | bothNaN);

    absDiff = abs(a - b);
    absDiff(bothNaN) = 0;
    scale = max(abs(a), abs(b));
    relDiff = absDiff ./ scale;
    relDiff(scale == 0) = 0;

    withinTolerance = differing & (absDiff <= options.AbsTol) & (relDiff <= options.RelTol);
    differing = differing & ~withinTolerance;

    entry.maxAbsDiff = max([0; absDiff(~bothNaN)]);
    entry.maxRelDiff = max([0; relDiff(~bothNaN)]);
    entry.nDiffering = sum(differing);
    entry.equal = entry.nDiffering == 0;
    if ~entry.equal
        entry.firstDiffIndex = find(differing, 1);
    end
else
    % Structs, cells, strings: compared whole rather than element-wise, so the
    % numeric difference columns stay empty rather than reporting a misleading
    % NaN. mvt.matHash still covers these contents, recursively.
    entry.equal = isequaln(expectedValue, actualValue);
    entry.note = sprintf('%s compared with isequaln', class(expectedValue));
end
end

% ---------------------------------------------------------------------------
function printReport(report, expectedFile, actualFile)
fprintf('\nComparing\n  expected: %s\n  actual:   %s\n\n', expectedFile, actualFile);

if ~isempty(report.onlyExpected)
    fprintf('  missing from actual: %s\n', strjoin(report.onlyExpected', ', '));
end
if ~isempty(report.onlyActual)
    fprintf('  unexpected in actual: %s\n', strjoin(report.onlyActual', ', '));
end

fprintf('%-34s %-8s %12s %12s %10s\n', 'VARIABLE', 'STATE', 'maxAbsDiff', 'maxRelDiff', 'nDiffer');
for iVar = 1:numel(report.variables)
    entry = report.variables(iVar);
    if entry.equal
        state = 'equal';
    else
        state = 'DIFFERS';
    end
    if isnan(entry.maxAbsDiff)
        fprintf('%-34s %-8s %12s %12s %10s', entry.name, state, '-', '-', '-');
    else
        fprintf('%-34s %-8s %12.4g %12.4g %10d', entry.name, state, ...
            entry.maxAbsDiff, entry.maxRelDiff, entry.nDiffering);
    end
    if ~isempty(entry.note)
        fprintf('  %s', entry.note);
    elseif ~entry.equal
        fprintf('  first at index %d', entry.firstDiffIndex);
    end
    fprintf('\n');
end

if report.identical
    fprintf('\nContents are identical.\n');
else
    fprintf('\nContents DIFFER.\n');
end
end

% ---------------------------------------------------------------------------
function options = parseOptions(varargin)
options = struct('AbsTol', 0, 'RelTol', 0, 'Verbose', true);
for iArg = 1:2:numel(varargin)
    name = validatestring(varargin{iArg}, {'AbsTol', 'RelTol', 'Verbose'});
    options.(name) = varargin{iArg + 1};
end
end
