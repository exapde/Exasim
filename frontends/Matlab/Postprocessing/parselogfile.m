function values = parselogfile(logfile, key)

%EXTRACT_TIMES  Extract numeric values following a given keyword from a log file.
%
%   values = extract_times('log.txt', 'hdgAssembleLinearSystem time:')
%
% Example:
%   extract_times('log.txt', 'Matrix-vector product time:')
%   extract_times('log.txt', 'Applying preconditioner time:')
%
% Parameters
% ----------
% logfile : string
%     Path to the log file.
% key : string
%     Text preceding the numeric value to extract.
%
% Returns
% -------
% values : double array
%     Vector of numeric values found after the specified key.

    key = char(key);

    % Read the full text from the log file
    txt = fileread(logfile);

    % Escape any regex special characters in key
    safeKey = regexptranslate('escape', key);

    % Build regex pattern: match the key followed by a number
    pattern = [safeKey, '\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+\-]?\d+)?)'];

    % Extract tokens
    tokens = regexp(txt, pattern, 'tokens');

    if isempty(tokens)
        warning('No entries found for key: "%s"', key);
        values = [];
        return;
    end
    
    % Convert cell array of tokens to numeric vector
    values = cellfun(@(c) str2double(c{1}), tokens(:));
end

