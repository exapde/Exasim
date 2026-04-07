function kernelTimes = read_ncu_kernel_times(csvFile)
%READ_NCU_KERNEL_TIMES  Read Nsight Compute CSV and aggregate GPU time per kernel.
%
% Usage:
%   T = read_ncu_kernel_times("ASMNCUraw.csv");
%   T = read_ncu_kernel_times("ASMNCUraw.cvs");   % (typo-safe; still tries to read)
%
% Output:
%   kernelTimes : table with per-kernel total time, mean time, and launch count,
%                 sorted by total time descending.

    if nargin < 1 || strlength(string(csvFile))==0
        csvFile = "ASMNCUraw.csv";
    end
    csvFile = string(csvFile);

    if ~isfile(csvFile)
        error("File not found: %s", csvFile);
    end

    % Robust CSV import (handles quoted fields, etc.)
    opts = detectImportOptions(csvFile, "Delimiter", ",");
    opts = setvaropts(opts, opts.VariableNames, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, opts.VariableNames, "EmptyFieldRule", "auto");
    raw = readtable(csvFile, opts);

    if isempty(raw) || width(raw)==0
        error("No data read from %s", csvFile);
    end

    % Normalize variable names for matching
    vnames = string(raw.Properties.VariableNames);
    vnorm  = lower(strrep(vnames, "_", " "));
    vnorm  = lower(strrep(vnorm, "-", " "));

    % -------- Find kernel name column --------
    % Common Nsight Compute exports use "Kernel Name" (sometimes with extra fields).
    nameCandidates = find( ...
        contains(vnorm, "kernel name") | ...
        (contains(vnorm, "kernel") & contains(vnorm, "name")) | ...
        contains(vnorm, "demangled name") | ...
        contains(vnorm, "function name") );

    if isempty(nameCandidates)
        error("Couldn't find a kernel name column. Columns are: %s", strjoin(vnames, ", "));
    end
    nameCol = nameCandidates(1);

    kernelName = raw.(vnames(nameCol));
    if ~isstring(kernelName), kernelName = string(kernelName); end
    kernelName = strtrim(kernelName);

    % -------- Find a GPU time column --------
    % Try to locate a per-kernel time metric column. Nsight CSVs vary a lot:
    % examples: "GPU Time", "Time", "Duration", "Kernel Time", "Avg", "Total Time", etc.
    timeCandidates = find( ...
        contains(vnorm, "gpu time") | ...
        contains(vnorm, "kernel time") | ...
        contains(vnorm, "duration") | ...
        (contains(vnorm, "time") & ~contains(vnorm, "timestamp")) );

    if isempty(timeCandidates)
        error("Couldn't find a time column. Columns are: %s", strjoin(vnames, ", "));
    end

    % Prefer the most specific-looking time column
    % (GPU Time > Kernel Time > Duration > Time)
    prio = zeros(size(timeCandidates));
    for k = 1:numel(timeCandidates)
        s = vnorm(timeCandidates(k));
        if contains(s,"gpu time"),        prio(k)=4;
        elseif contains(s,"kernel time"), prio(k)=3;
        elseif contains(s,"duration"),    prio(k)=2;
        else,                            prio(k)=1;
        end
    end
    [~,ix] = max(prio);
    timeCol = timeCandidates(ix);

    timeRaw = raw.(vnames(timeCol));

    % Convert time column to numeric seconds (auto-detect units if embedded as text)
    timeSec = convert_time_to_seconds(timeRaw);

    % Drop rows with empty kernel names or NaN times
    ok = kernelName ~= "" & ~isnan(timeSec);
    kernelName = kernelName(ok);
    timeSec    = timeSec(ok);

    % Aggregate
    [G, kernel] = findgroups(kernelName);
    totalSec    = splitapply(@sum,  timeSec, G);
    meanSec     = splitapply(@mean, timeSec, G);
    count       = splitapply(@numel,timeSec, G);

    kernelTimes = table(kernel, count, totalSec, meanSec, ...
        'VariableNames', {'Kernel','Launches','TotalTime_s','MeanTime_s'});
    kernelTimes = sortrows(kernelTimes, 'TotalTime_s', 'descend');

    % Optional: also show ms for convenience
    kernelTimes.TotalTime_ms = 1e3 * kernelTimes.TotalTime_s;
    kernelTimes.MeanTime_ms  = 1e3 * kernelTimes.MeanTime_s;

    % Print top 20
    disp(kernelTimes(1:min(20,height(kernelTimes)), :));

    % Write summary next to input
    [p,n,~] = fileparts(csvFile);
    outFile = fullfile(p, n + "_kernel_times.csv");
    writetable(kernelTimes, outFile);
    fprintf("Wrote summary: %s\n", outFile);
end

function tsec = convert_time_to_seconds(x)
%CONVERT_TIME_TO_SECONDS Convert numeric or string time values to seconds.
% Handles cases like:
%   123.4              (assumed already seconds)
%   "12.3 ms", "7 us", "0.5 s", "100 ns"
%   "1,234.5"          (commas)
% If your CSV is known to be milliseconds, you can hardcode conversion.

    if isnumeric(x)
        t = double(x);
        % Heuristic: if typical values look like milliseconds (e.g., > 1e3), you can adapt.
        tsec = t; % assume seconds by default
        return;
    end

    xs = string(x);
    xs = strtrim(xs);
    xs = replace(xs, ",", "");   % remove thousands separators

    % Extract numeric part and unit (if any)
    num = nan(size(xs));
    unit = strings(size(xs));

    for i = 1:numel(xs)
        s = xs(i);
        if s == "" || s == "NaN"
            continue;
        end
        % Match: number + optional unit
        tok = regexp(s, "^\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([a-zA-Zµ]*)\s*$", "tokens", "once");
        if isempty(tok)
            continue;
        end
        num(i)  = str2double(tok{1});
        unit(i) = lower(string(tok{2}));
        if unit(i) == ""
            unit(i) = "s"; % assume seconds if unit missing
        end
        unit(i) = replace(unit(i), "us", "µs"); % normalize
    end

    tsec = nan(size(num));
    for i = 1:numel(num)
        if isnan(num(i)), continue; end
        switch unit(i)
            case {"s","sec","secs","second","seconds"}
                tsec(i) = num(i);
            case {"ms","msec","msecs","millisecond","milliseconds"}
                tsec(i) = num(i) * 1e-3;
            case {"us","µs","usec","usecs","microsecond","microseconds"}
                tsec(i) = num(i) * 1e-6;
            case {"ns","nsec","nsecs","nanosecond","nanoseconds"}
                tsec(i) = num(i) * 1e-9;
            otherwise
                % Unknown unit -> assume seconds
                tsec(i) = num(i);
        end
    end
end