function digest = matHash(source, varargin)
% MVT.MATHASH  Content hash of a .mat file's variables, ignoring container noise.
%
% Purpose
%   Byte hashes work for the pipeline's JSON outputs but not for its .mat
%   outputs: v7 files are gzip streams carrying a creation timestamp, and v7.3
%   files are HDF5 containers whose layout depends on the writer. Two saves of
%   identical data therefore differ on disk. This hashes the *contents* -
%   variable names, classes, sizes, and values - in a canonical order, so the
%   result depends only on the data.
%
% Inputs
%   source    path to a .mat file, or a struct of already-loaded variables
%   varargin  name/value pairs:
%               'Decimals'  round numeric data to this many decimals before
%                           hashing; [] (default) hashes exact bit patterns
%               'Exclude'   cellstr of variable names to skip
%
% Outputs
%   digest  char, lowercase hex MD5 of the canonical serialization
%
% Algorithm
%   1. Load the file (or accept the struct).
%   2. Visit variables in sorted name order; for each, feed the digest the
%      name, class, and size, then the data:
%        - numeric/logical/char: raw bytes via typecast (optionally rounded)
%        - struct: field names in sorted order, then each field, recursively
%        - cell: elements in column-major order, recursively
%        - anything else: its char representation
%   3. Return the digest as hex.
%
% Notes
%   Rounding is the knob to use when comparing across MATLAB releases, where
%   the last bits of transcendental functions may legitimately differ:
%     mvt.matHash(file, 'Decimals', 6)
%
% Dependencies
%   None beyond MATLAB (uses the bundled Java MessageDigest).
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

options = parseOptions(varargin{:});

if ischar(source) || isstring(source)
    source = char(source);
    if ~isfile(source)
        error('mvt:matHash:missingFile', 'No such file: %s', source);
    end
    data = load(source);
else
    data = source;
end

names = sort(fieldnames(data));
names = names(~ismember(names, options.Exclude));

engine = java.security.MessageDigest.getInstance('MD5');
for iName = 1:numel(names)
    updateText(engine, sprintf('var:%s', names{iName}));
    hashValue(engine, data.(names{iName}), options);
end

digest = lower(reshape(dec2hex(typecast(engine.digest(), 'uint8')).', 1, []));
end

% ---------------------------------------------------------------------------
function hashValue(engine, value, options)
updateText(engine, sprintf('class:%s|size:%s', class(value), mat2str(size(value))));

if isnumeric(value) || islogical(value)
    if isempty(value)
        return
    end
    if isa(value, 'double') && ~isempty(options.Decimals)
        value = round(value, options.Decimals);
    end
    if islogical(value)
        value = uint8(value);
    end
    engine.update(typecast(value(:), 'uint8'));

elseif ischar(value)
    updateText(engine, value(:).');

elseif isstruct(value)
    fields = sort(fieldnames(value));
    for iElement = 1:numel(value)
        for iField = 1:numel(fields)
            updateText(engine, sprintf('field:%s', fields{iField}));
            hashValue(engine, value(iElement).(fields{iField}), options);
        end
    end

elseif iscell(value)
    for iElement = 1:numel(value)
        hashValue(engine, value{iElement}, options);
    end

else
    % Objects (datetime, categorical, ...) hash through their text form.
    updateText(engine, char(string(value(:).')));
end
end

% ---------------------------------------------------------------------------
function updateText(engine, text)
engine.update(uint8(char(text)));
end

% ---------------------------------------------------------------------------
function options = parseOptions(varargin)
options = struct('Decimals', [], 'Exclude', {{}});
for iArg = 1:2:numel(varargin)
    name = validatestring(varargin{iArg}, {'Decimals', 'Exclude'});
    options.(name) = varargin{iArg + 1};
end
if ischar(options.Exclude) || isstring(options.Exclude)
    options.Exclude = cellstr(options.Exclude);
end
end
