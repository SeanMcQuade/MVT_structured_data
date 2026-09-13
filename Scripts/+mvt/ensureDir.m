function ensureDir(varargin)
% MVT.ENSUREDIR  Create directories if they do not already exist.
%
% Purpose
%   Replaces the repeated `if ~isfolder(p), mkdir(p), end` blocks in the stage
%   scripts, and makes directory creation safe when several sharded MATLAB
%   processes start at the same moment (mkdir can lose a race, so a second
%   isfolder check decides success).
%
% Inputs
%   One or more folder paths (char, string, or cellstr).
%
% Outputs
%   (none; raises mvt:ensureDir:failed if a folder cannot be created)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

for iArg = 1:numel(varargin)
    arg = varargin{iArg};
    if isstring(arg) || ischar(arg)
        arg = cellstr(arg);
    end
    for iPath = 1:numel(arg)
        folder = char(arg{iPath});
        if isempty(folder) || isfolder(folder)
            continue
        end
        [ok, msg] = mkdir(folder);
        if ~ok && ~isfolder(folder)
            % isfolder() re-checked: a concurrent worker may have won the race.
            error('mvt:ensureDir:failed', ...
                'Could not create folder %s: %s', folder, msg);
        end
    end
end
end
