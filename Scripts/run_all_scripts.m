function [] = run_all_scripts(varargin)
% RUN_ALL_SCRIPTS  Run the whole MVT pipeline from inside MATLAB.
%
% Purpose
%   Generates the slim and full versions of the MVT data for the three test
%   days (16-18/11/2022), aggregates samples from the data, builds the
%   macroscopic fields, and generates the results figures - the same graph of
%   work the Makefile drives, for people who would rather stay in MATLAB.
%
% Inputs
%   varargin  options struct and/or name/value pairs (see mvt.options), e.g.
%               run_all_scripts                       % all days, skip fresh work
%               run_all_scripts('Days', 17)           % one day
%               run_all_scripts('Force', true)        % rebuild everything
%               run_all_scripts('DryRun', true)       % show the plan only
%
% Outputs
%   (none; each stage writes its own files under results/)
%
% Notes
%   Every stage checks for itself whether its outputs are missing, older than
%   their inputs, or older than the code that produced them, so re-running this
%   is cheap. For parallel execution across days and segments use the Makefile:
%     make -j3 all           # three days at once
%     make SHARDS=4 slim-17  # four MATLAB processes over one day's segments
%
% Dependencies
%   mvt.options, mvt.build
%
% (C) 2025-2026 CIRCLES Consortium. Author: Sulaiman Almatrudi. BSD-3-Clause.

opts = mvt.options(varargin{:});

% Run through each of the test days
for processingDay = opts.Days
    % Assemble GPS data into run trajectories
    mvt.build('gps', processingDay, opts);
    % Generate full version of the data for a full test day
    mvt.build('full', processingDay, opts);
    % Generate slim version of the data for a full test day
    mvt.build('slim', processingDay, opts);
    % Generate samples for analysis
    mvt.build('samples', processingDay, opts);
    % Generate macroscopic fields from MOTION data
    mvt.build('fields', processingDay, opts);
    % Plot macroscopic fields
    mvt.build('macro', processingDay, opts);
    % Plot microscopic trajectories
    mvt.build('micro', processingDay, opts);
end

% Plot stats of controlled AVs and sampled leaders/followers, comparing days
mvt.build('av', [], opts);
end
