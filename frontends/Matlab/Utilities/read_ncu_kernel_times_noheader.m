function kernelTimes = read_ncu_kernel_times_noheader(csvFile, kernelCol, timeCol)
%READ_NCU_KERNEL_TIMES_NOHEADER  Aggregate GPU time per kernel from Nsight CSV,
%even when MATLAB imports columns as Var1..VarN (no recognized header).
%
% Usage:
%   T = read_ncu_kernel_times_noheader("ASMNCUraw.csv");
%   T = read_ncu_kernel_times_noheader("ASMNCUraw.csv", 5, 12); % manual override
%
% Optional:
%   kernelCol : column index containing kernel names
%   timeCol   : column index containing per-launch time (any unit; see notes)
%
% Output:
%   kernelTimes : table sorted by total time (descending)

    if nargin < 1 || strlength(string(csvFile))==0
        csvFile = "ASMNCUraw.csv";
    end
    csvFile = string(csvFile);

    if ~isfile(csvFile)
        error("File not found: %s", csvFile);
    end

    % Read EVERYTHING as raw cells (robust to weird headers/preamble)
    C = readcell(csvFile, "Delimiter", ",");

    if isempty(C) || size(C,2) < 2
        error("Could not read usable data from %s", csvFile);
    end

    % Drop completely empty rows
    emptyRow = false(size(C,1),1);
    for i=1:size(C,1)
        row = C(i,:);
        emptyRow(i) = all(cellfun(@(v) (isempty(v) || (isstring(v)&&strlength(v)==0) || (ischar(v)&&isempty(strtrim(v)))), row));
    end
    C(emptyRow,:) = [];

    % Convert to string for type checks
    S = cellfun(@cell_to_string, C, "UniformOutput", false);
    S = string(vertcat(S{:}));
    S = reshape(S, size(C,1), size(C,2));

    % Find first row that looks like "data" (has some numeric columns)
    dataStart = 1;
    for r = 1:size(C,1)
        numCount = 0;
        for c = 1:size(C,2)
            if is_numeric_like(S(r,c))
                numCount = numCount + 1;
            end
        end
        if numCount >= 2
            dataStart = r;
            break;
        end
    end

    D = C(dataStart:end, :);
    DS = S(dataStart:end, :);

    nrow = size(D,1);
    ncol = size(D,2);

    % ---- Auto-detect kernel and time columns if not provided ----
    if nargin < 2 || isempty(kernelCol) || nargin < 3 || isempty(timeCol)
        [kGuess, tGuess] = guess_kernel_and_time_columns(DS);
        if nargin < 2 || isempty(kernelCol), kernelCol = kGuess; end
        if nargin < 3 || isempty(timeCol),   timeCol   = tGuess; end
    end

    if kernelCol < 1 || kernelCol > ncol
        error("kernelCol=%d out of range (1..%d).", kernelCol, ncol);
    end
    if timeCol < 1 || timeCol > ncol
        error("timeCol=%d out of range (1..%d).", timeCol, ncol);
    end

    kernelName = strtrim(DS(:, kernelCol));
    tsec = convert_time_to_seconds(DS(:, timeCol));

    ok = kernelName ~= "" & ~isnan(tsec);
    kernelName = kernelName(ok);
    tsec       = tsec(ok);

    if isempty(kernelName)
        error("No kernel names parsed. Try specifying kernelCol manually.");
    end
    if all(isnan(tsec))
        error("No time values parsed. Try specifying timeCol manually.");
    end

    [G, kernel] = findgroups(kernelName);
    totalSec = splitapply(@sum,  tsec, G);
    meanSec  = splitapply(@mean, tsec, G);
    count    = splitapply(@numel,tsec, G);

    kernelTimes = table(kernel, count, totalSec, meanSec, ...
        'VariableNames', {'Kernel','Launches','TotalTime_s','MeanTime_s'});
    kernelTimes.TotalTime_ms = 1e3 * kernelTimes.TotalTime_s;
    kernelTimes.MeanTime_ms  = 1e3 * kernelTimes.MeanTime_s;

    kernelTimes = sortrows(kernelTimes, 'TotalTime_s', 'descend');

    disp(kernelTimes(1:min(20,height(kernelTimes)), :));

    [p,n,~] = fileparts(csvFile);
    outFile = fullfile(p, n + "_kernel_times.csv");
    writetable(kernelTimes, outFile);
    fprintf("Wrote summary: %s\n", outFile);
    fprintf("Used kernelCol=%d, timeCol=%d (data started at row %d of original file)\n", ...
        kernelCol, timeCol, dataStart);
end

% ---------------- helpers ----------------

function s = cell_to_string(v)
    if isstring(v), s = v;
    elseif ischar(v), s = string(v);
    elseif isnumeric(v) && isscalar(v), s = string(v);
    elseif ismissing(v), s = "";
    else, s = string(v);
    end
end

function tf = is_numeric_like(s)
    s = strtrim(s);
    if s=="" || lower(s)=="nan"
        tf = false; return;
    end
    s2 = replace(s, ",", "");
    tf = ~isnan(str2double(s2));
end

function [kernelCol, timeCol] = guess_kernel_and_time_columns(DS)
    % Heuristic:
    %  - kernel column: mostly non-numeric strings, high uniqueness, medium/long tokens
    %  - time column: mostly numeric-like or numeric-with-units, low-ish magnitude
    ncol = size(DS,2);
    nrow = size(DS,1);

    nonNumFrac = zeros(1,ncol);
    uniqFrac   = zeros(1,ncol);
    avgLen     = zeros(1,ncol);
    hasAngleOrColon = zeros(1,ncol);

    for c = 1:ncol
        col = DS(:,c);
        isNum = arrayfun(@(x) is_numeric_like(x), col);
        nonNumFrac(c) = 1 - mean(isNum);

        col2 = col(col~="");
        uniqFrac(c) = numel(unique(col2)) / max(1,numel(col2));
        avgLen(c) = mean(strlength(col2));

        hasAngleOrColon(c) = mean(contains(col2,"::") | contains(col2,"<") | contains(col2,">") | contains(col2,"("));
    end

    % Kernel score: want high non-numeric fraction + decent uniqueness + longer strings
    kernelScore = 2.0*nonNumFrac + 1.0*uniqFrac + 0.02*avgLen + 0.5*hasAngleOrColon;
    [~, kernelCol] = max(kernelScore);

    % Time column guess: among columns that look numeric-ish, choose one with
    % high numeric fraction and non-trivial variance
    numFrac = 1 - nonNumFrac;
    timeScore = -inf(1,ncol);

    for c = 1:ncol
        if numFrac(c) < 0.6, continue; end
        t = convert_time_to_seconds(DS(:,c));
        t = t(~isnan(t));
        if numel(t) < max(10, 0.1*nrow), continue; end
        timeScore(c) = numFrac(c) + 0.1*log(1+std(t));
    end

    [best, timeCol] = max(timeScore);
    if ~isfinite(best)
        % fallback: pick most numeric-looking column (excluding kernelCol)
        numFrac(kernelCol) = -inf;
        [~, timeCol] = max(numFrac);
    end
end

function tsec = convert_time_to_seconds(col)
    % Accepts string array; supports:
    %   "12.3 ms", "7 us", "0.5 s", "100 ns", or plain numbers as strings.
    xs = string(col);
    xs = strtrim(xs);
    xs = replace(xs, ",", "");

    tsec = nan(size(xs));
    for i = 1:numel(xs)
        s = xs(i);
        if s=="" || lower(s)=="nan"
            continue;
        end
        tok = regexp(s, "^\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([a-zA-Zµ]*)\s*$", "tokens", "once");
        if isempty(tok), continue; end

        val = str2double(tok{1});
        unit = lower(string(tok{2}));
        if unit=="", unit="s"; end
        unit = replace(unit, "us", "µs");

        switch unit
            case {"s","sec","secs","second","seconds"}
                tsec(i) = val;
            case {"ms","msec","msecs","millisecond","milliseconds"}
                tsec(i) = val * 1e-3;
            case {"us","µs","usec","usecs","microsecond","microseconds"}
                tsec(i) = val * 1e-6;
            case {"ns","nsec","nsecs","nanosecond","nanoseconds"}
                tsec(i) = val * 1e-9;
            otherwise
                % Unknown unit -> assume seconds
                tsec(i) = val;
        end
    end
end