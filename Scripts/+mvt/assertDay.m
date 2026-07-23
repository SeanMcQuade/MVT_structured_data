function assertDay(day)
% MVT.ASSERTDAY  Validate a MegaVanderTest processing day.
%
% Purpose
%   The experiment ran on 16, 17, and 18 November 2022. Every stage takes the
%   day as a bare number; this check turns a typo into an immediate, readable
%   error instead of an empty directory listing several minutes later.
%
% Inputs
%   day  numeric scalar; must be 16, 17, or 18
%
% Outputs
%   (none; raises mvt:assertDay:badDay on failure)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if ~isnumeric(day) || ~isscalar(day) || ~ismember(day, [16 17 18])
    error('mvt:assertDay:badDay', ...
        ['Specify the day of Nov. 2022 MVT to process (16, 17, or 18); ', ...
        'got %s.'], mat2str(day));
end
end
