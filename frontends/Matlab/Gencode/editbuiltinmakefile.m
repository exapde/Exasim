function editbuiltinmakefile(modelid, makefile)
%EDITBUILTINMAKEFILE Append model<id>/model.cpp to SRC line if missing
%
%   editbuiltinmakefile(3)
%   editbuiltinmakefile(3, 'Makefile')

    if nargin < 2
        makefile = 'Makefile';
    end

    if ~isfile(makefile)
        error('Makefile not found: %s', makefile);
    end

    % Read entire Makefile as text
    txt = fileread(makefile);

    % Construct the model source string
    newModel = sprintf('model%d/model.cpp', modelid);

    % If already present, do nothing
    if ~isempty(strfind(txt, newModel)) %#ok<STREMP>
        fprintf('model%d already present in SRC. No change made.\n', modelid);
        return;
    end

    % Find the SRC line
    key = 'SRC :=';
    k = strfind(txt, key);
    if isempty(k)
        error('SRC := line not found in Makefile');
    end

    % Extract the full SRC line
    lineStart = k(1);
    lineEnd = strfind(txt(lineStart:end), newline);
    lineEnd = lineStart + lineEnd(1) - 2;

    oldLine = txt(lineStart:lineEnd);

    % Append new model source
    newLine = [oldLine ' ' newModel];

    % Replace in file
    txt = strrep(txt, oldLine, newLine);

    % Write back to file
    fid = fopen(makefile, 'w');
    fwrite(fid, txt);
    fclose(fid);

    fprintf('Added %s to SRC in %s\n', newModel, makefile);
end