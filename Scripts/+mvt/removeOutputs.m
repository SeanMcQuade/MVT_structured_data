function removed = removeOutputs(outputs, opts)
% MVT.REMOVEOUTPUTS  Delete a stage's outputs (the Clean option / make clean).
%
% Purpose
%   Implements `Clean` for stages whose outputs are described by globs, such as
%   the figure scripts. Deliberately narrow: it deletes only files that match
%   the patterns a stage declares in mvt.expectedOutputs, never a whole folder,
%   because results/ holds 150 GB that must not be removed by accident.
%
% Inputs
%   outputs  char / cellstr of paths or globs (from mvt.expectedOutputs)
%   opts     (optional) options struct; DryRun reports without deleting,
%            Verbose logs each deletion
%
% Outputs
%   removed  cellstr of the files that were (or would be) deleted
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if nargin < 2 || isempty(opts)
    opts = mvt.options();
end
if ischar(outputs) || isstring(outputs)
    outputs = cellstr(outputs);
end

removed = {};
for iOut = 1:numel(outputs)
    spec = char(outputs{iOut});
    folder = fileparts(spec);
    listing = dir(spec);
    listing = listing(~[listing.isdir]);
    for iFile = 1:numel(listing)
        target = fullfile(folder, listing(iFile).name);
        if opts.DryRun
            mvt.log(opts, 'would delete %s', target);
        else
            delete(target);
            mvt.log(opts, 'deleted %s', target);
        end
        removed{end+1} = target; %#ok<AGROW>
    end
end
end
