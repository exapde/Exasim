function strArray = convertHandlesToStrings(fhandleArray)
% Convert a cell array of function handles into descriptive strings

    n = numel(fhandleArray);
    strArray = strings(1, n);
    
    for i = 1:n
        fstr = func2str(fhandleArray{i});  % e.g., '@(p)abs(p(2,:))<1e-8'
        
        % Remove '@(p)' prefix
        fstr = erase(fstr, '@(p)');
        
        % Replace p(1,:) with x, p(2,:) with y
        fstr = regexprep(fstr, 'p\(1,:\)', 'x');
        fstr = regexprep(fstr, 'p\(2,:\)', 'y');
        fstr = regexprep(fstr, 'p\(3,:\)', 'z');
        
        fstr = replace(fstr,"p([1,2],:)", "xy");
        fstr = replace(fstr,"p([1,3],:)", "xz");
        fstr = replace(fstr,"p([2,3],:)", "yz");
        
        % Optional: clean up formatting
        fstr = strtrim(fstr);
        
        strArray(i) = fstr;
    end
end



