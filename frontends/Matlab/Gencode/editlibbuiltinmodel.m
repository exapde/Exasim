function editlibbuiltinmodel(modelid, cppfile)
%EDITLIBBUILTINMODEL Ensure libbuiltinmodel.cpp dispatches builtin model "modelid".
%
%   editlibbuiltinmodel(modelid)
%   editlibbuiltinmodel(modelid, cppfile)
%
% Idempotent: re-running will not duplicate includes/cases.
    
    if ~isfile(cppfile)
        error("File not found: %s", cppfile);
    end

    txt = readTextFile(cppfile);
    
    mid  = modelid;
    inc  = sprintf('#include "model%d/model.hpp"\n', mid);
    nsln = sprintf('namespace m%d = exasim_model_%d;\n', mid, mid);
                
    % 1) Ensure include line exists
    if ~contains(txt, inc) && ~contains(txt, nsln)
        txt = insertModelInclude(txt, inc);
        txt = insertModelNamespace(txt, nsln);
        txt = insertCasesEverywhere(txt, mid);
        writeTextFile(cppfile, txt);        
        fprintf("Updated %s to include model%d.\n", cppfile, mid);
    end
    
end

% ---------------- helpers ----------------

function txt = readTextFile(fname)
    txt = fileread(fname);
end

function writeTextFile(fname, txt)
    fid = fopen(fname, 'w');
    if fid < 0
        error("Cannot open for writing: %s", fname);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, txt, 'char');
end

function txt = insertModelInclude(txt, incLine)
    % Insert after the last: #include "model<number>/model.hpp"
    pat = '#include\s+"model\d+/model\.hpp"\s*';
    [s,e] = regexp(txt, pat, 'start', 'end');

    if ~isempty(e)
        insertPos = e(end);
        txt = txtInsert(txt, insertPos, incLine);
        return;
    end

    % Fallback: insert after "libbuiltinmodel.hpp" include if present
    pat2 = '#include\s+"libbuiltinmodel\.hpp"\s*';
    [s2,e2] = regexp(txt, pat2, 'start', 'end', 'once');
    if ~isempty(e2)
        insertPos = e2;
        txt = txtInsert(txt, insertPos, incLine );
        return;
    end

    % Last resort: after the last #include
    pat3 = '#include[^\n]*\n';
    [~,e3] = regexp(txt, pat3, 'start', 'end');
    if ~isempty(e3)
        insertPos = e3(end);
        txt = txtInsert(txt, insertPos, incLine);
    else
        % If no includes, prepend
        txt = incLine + txt;
    end
end

function txt = insertModelNamespace(txt, nsLine)
    % Insert after the last namespace alias line:
    % namespace mX = exasim_model_X;
    pat = 'namespace\s+m\d+\s*=\s*exasim_model_\d+;\s*';
    [~,e] = regexp(txt, pat, 'start', 'end');

    if ~isempty(e)
        insertPos = e(end);
        txt = txtInsert(txt, insertPos, nsLine);
        return;
    end

    % Fallback: insert after the last include block
    pat2 = '#include[^\n]*\n';
    [~,e2] = regexp(txt, pat2, 'start', 'end');
    if ~isempty(e2)
        insertPos = e2(end);
        txt = txtInsert(txt, insertPos, nsLine);
    else
        % Prepend if needed
        txt = nsLine + txt;
    end
end

function txt = insertCasesEverywhere(txt, mid)
    % For each "switch (builtinmodelID) { ... default:" block, insert
    % "case mid: m<mid>::<FuncName>(...); return;" before default if missing.

    % Find all function blocks that look like:
    % extern "C" \n void NAME( ... ) \n { \n switch (builtinmodelID) { ... } \n }
    % We'll do a simpler approach:
    % - Locate each "extern "C"\nvoid <fname>("
    % - Within that function body, find "switch (builtinmodelID)" and "default:"
    % - Determine function name (<fname>)
    %
    % This avoids trying to parse the entire C++ grammar.

    funcPat = 'extern\s+"C"\s*\n\s*void\s+(\w+)\s*\(';
    [starts, ends, tokens] = regexp(txt, funcPat, 'start', 'end', 'tokens');

    if isempty(starts)
        warning("No extern ""C"" void functions found; no switch blocks updated.");
        return;
    end
    
    % Process from end to start to keep indices stable when inserting
    for k = numel(starts):-1:1
        fname = tokens{k}{1};

        % Determine a conservative slice of the function body:
        % from the end of the signature match to the next 'extern "C"' or end-of-file
        bodyStart = ends(k) + 1;
        if k < numel(starts)
            bodyEnd = starts(k+1) - 1;
        else
            bodyEnd = strlength(txt);
        end

        block = extractBetween(txt, bodyStart, bodyEnd);
        
        % Only touch functions that have switch(builtinmodelID)
        if isempty(regexp(block, 'switch\s*\(\s*builtinmodelID\s*\)\s*\{', 'once'))
            continue;
        end

        % If case already exists, skip
        casePat = sprintf('case\\s+%d\\s*:\\s*m%d::%s\\s*\\(', mid, mid, fname);
        if ~isempty(regexp(block, casePat, 'once'))
            continue;
        end

        % Find "default:" inside this function block (first occurrence after switch)
        defIdx = regexp(block, '\n\s*default\s*:', 'once');
        if isempty(defIdx)
            warning("No default: found in function %s; skipping.", fname);
            continue;
        end

        % Build insertion line matching your one-line style
        insertLine = sprintf('%s      case %d: m%d::%s(%s); return;', ...
            newline, mid, mid, fname, '/* args */');

        % Replace /* args */ with the argument list used in an existing case
        insertLine = fillArgsFromExistingCase(block, fname, insertLine);
        
        % Insert right before default:
        blockNew = block;
        blockNew = insertBefore(blockNew, defIdx, insertLine);

        % Put back
        txt = replaceBetween(txt, bodyStart, bodyEnd, blockNew);
    end
end

function insertLine = fillArgsFromExistingCase(block, fname, insertLine)
    % Extract the argument list from any existing "case 1:" line calling this function,
    % so we can replicate exactly.
    %
    % Example line:
    % case 1: m1::KokkosFlux(f, xdg, ...); return;
    callPat = sprintf('case\\s+\\d+\\s*:\\s*m\\d+::%s\\s*\\(([^;]*)\\)\\s*;\\s*return\\s*;', fname);
    tok = regexp(block, callPat, 'tokens', 'once');

    if isempty(tok)
        % If not found, fall back to using the function signature in the file is hard;
        % keep placeholder so user notices. Better than silently wrong.
        return;
    end

    args = strtrim(tok{1}); % inside parentheses
    insertLine = strrep(insertLine, '/* args */', args);
end

function out = txtInsert(txt, pos, insertStr)
    % pos is 1-based insertion index into the character vector/string
    txt = char(txt); % work as char for indexing
    if pos < 1, pos = 1; end
    if pos > numel(txt)+1, pos = numel(txt)+1; end
    out = [txt(1:pos-1), char(insertStr), txt(pos:end)];
    out = string(out);
end