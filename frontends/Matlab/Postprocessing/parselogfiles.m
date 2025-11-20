function [data, logfiles] = parselogfiles(rootdir, ext)
%parselogfiles parse all log files (e.g. *.out, *.log) recursively.
%
%   txt = parselogfiles('/path/to/folder', '*.out')
%   txt = parselogfiles('C:\Users\me\logs', '*.log')

    if nargin < 2
        ext = '*.out'; % default
    end

    % Find all matching files recursively
    files = dir(fullfile(rootdir, '**', ext));

    if isempty(files)
        error('No files matching "%s" found under %s', ext, rootdir);
    end

    data = cell(numel(files), 1);
    logfiles = cell(numel(files), 1);
    for k = 1:numel(files)
        fname = fullfile(files(k).folder, files(k).name);
        fprintf('parsing %s\n', fname);
        data{k} = parselogfile1(fname);
        logfiles{k} = fname;
    end
end
