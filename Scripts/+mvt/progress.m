function report = progress(total, label, varargin)
% MVT.PROGRESS  Progress reporter for a stage's inner loop.
%
% Purpose
%   Stage output used to be a wall of "Loading ... Done (12sec)." lines with no
%   sense of how far along a run was. This returns a closure you call once per
%   item; it prints a bar, the item count, elapsed time, and an ETA derived from
%   the items completed so far.
%
% Usage
%   report = mvt.progress(numel(files), 'slim 2022-11-17');
%   for k = 1:numel(files)
%       ... do the work ...
%       report(k, files(k).name);       % or report(k) with no detail
%   end
%   report();                            % finish: prints the summary line
%
% Inputs
%   total  number of items
%   label  char, what is being processed (stage and day)
%   Name/value:
%     'Opts'    options struct; Verbose=false silences the reporter entirely
%     'Stream'  fid to write to (default 1, stdout)
%
% Outputs
%   report  function handle, report(k, detail) / report(k) / report()
%
% Notes
%   `matlab -batch` writes to a pipe under make, not a terminal, so this never
%   uses carriage returns or escape codes: every line is durable and greppable,
%   and interleaved output from `make -j` stays readable because each line
%   carries its own label. Lines are emitted at most once per UpdateSeconds so a
%   24-segment stage produces a couple of dozen lines, not thousands.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

parser = inputParser;
parser.addParameter('Opts', mvt.options());
parser.addParameter('Stream', 1);
parser.addParameter('UpdateSeconds', 5);
parser.parse(varargin{:});
opts = parser.Results.Opts;
stream = parser.Results.Stream;
updateEvery = parser.Results.UpdateSeconds;

quiet = ~(isstruct(opts) && isfield(opts, 'Verbose') && opts.Verbose);
started = tic;
lastPrint = -inf;
tag = shardTag(opts);

report = @step;

    function step(k, detail)
        if quiet
            return
        end
        if nargin == 0                       % final summary
            fprintf(stream, '%s %s | %s %d/%d done in %s\n', tag, label, ...
                bar(1, 20), total, total, humanTime(toc(started)));
            return
        end
        elapsed = toc(started);
        isLast = (k >= total);
        if ~isLast && (elapsed - lastPrint) < updateEvery
            return                            % rate-limit; keep logs readable
        end
        lastPrint = elapsed;
        fraction = k / max(total, 1);
        line = sprintf('%s %s | %s %d/%d | %s elapsed', tag, label, ...
            bar(fraction, 20), k, total, humanTime(elapsed));
        if k > 0 && ~isLast
            remaining = elapsed / k * (total - k);
            line = [line sprintf(' | ~%s left', humanTime(remaining))];
        end
        if nargin > 1 && ~isempty(detail)
            line = [line ' | ' char(detail)];
        end
        fprintf(stream, '%s\n', line);
    end
end

% ---------------------------------------------------------------------------
function text = bar(fraction, width)
fraction = min(1, max(0, fraction));
filled = round(fraction * width);
text = ['[' repmat('#', 1, filled) repmat('.', 1, width - filled) ']'];
end

% ---------------------------------------------------------------------------
function text = humanTime(seconds)
seconds = max(0, seconds);
if seconds < 60
    text = sprintf('%ds', round(seconds));
elseif seconds < 3600
    text = sprintf('%dm%02ds', floor(seconds/60), round(mod(seconds,60)));
else
    text = sprintf('%dh%02dm', floor(seconds/3600), floor(mod(seconds,3600)/60));
end
end

% ---------------------------------------------------------------------------
function tag = shardTag(opts)
tag = '[mvt]';
if isstruct(opts) && isfield(opts, 'Shard') && numel(opts.Shard) == 2 ...
        && opts.Shard(2) > 1
    tag = sprintf('[mvt %d/%d]', opts.Shard(1), opts.Shard(2));
end
end
