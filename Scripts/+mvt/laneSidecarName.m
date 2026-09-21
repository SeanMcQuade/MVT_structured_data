function name = laneSidecarName(segmentOutputName)
% MVT.LANESIDECARNAME  Sidecar filename pairing with one slim segment.
%
% Purpose
%   generate_orig_dist_lanes writes one origin/destination-lane sidecar per raw
%   segment, and extract_lane_changes_v_dist_to_av reads them back paired with
%   the slim JSON of the same segment. Deriving that name in one place keeps
%   the producer, the consumer, and mvt.expectedOutputs from drifting apart.
%
% Inputs
%   segmentOutputName  char, the processed segment name from mvt.manifest /
%                      mvt.segmentName, e.g. 'I-24MOTION_2022-11-16_05-59-59.json'
%
% Outputs
%   name  char, e.g. 'I-24MOTION_2022-11-16_05-59-59_orig_dist_lane.mat'
%
% Dependencies
%   (none)
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

segmentOutputName = char(segmentOutputName);
[~, base, ext] = fileparts(segmentOutputName);
if ~strcmpi(ext, '.json')
    error('mvt:laneSidecarName:notASegment', ...
        'Expected a segment .json name, got ''%s''.', segmentOutputName);
end
name = [base '_orig_dist_lane.mat'];
end
