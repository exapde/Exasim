function [files_sorted,order] = reorder_files(files)
%REORDER_FILES Reorder file paths by type, polynomial order, mesh size, and variant (0 or 10)
%
%   files_sorted = reorder_files(files)
%
%   Sorting rules:
%     1. "BJ" entries come before "ASM"
%     2. Within each type, sort by polynomial order: P1, P2, P3, P4
%     3. Within each polynomial order, sort by mesh size: 16, 32, 64, 128, 256
%     4. Within each mesh size, sort "BJ0" < "BJ10" and "ASM0" < "ASM10"

    % Convert nested cell array if needed
    if iscell(files{1})
        files = vertcat(files{:});
    end

    % Extract polynomial order (P1, P2, P3, P4)
    polyOrderStr = regexp(files, '/P(\d+)/', 'tokens', 'once');
    polyOrder = cellfun(@(x) str2double(x{1}), polyOrderStr);

    % % Extract mesh size (16, 32, 64, 128, 256)
    meshSizeStr = regexp(files, '/(\d{2,3})/', 'tokens', 'once');    
    meshsz = 1;
    if numel(meshSizeStr{1})==0
        meshsz = 0; 
    end

    % Determine category: 1 = BJ, 2 = ASM
    category = zeros(size(files));
    category(contains(files, 'BJ')) = 1;
    category(contains(files, 'ASM')) = 2;

    % Determine sub-type: 0 for *0, 10 for *10
    subtypeStr = regexp(files, '(BJ|ASM)(\d+)', 'tokens', 'once');
    subtype = cellfun(@(x) str2double(x{2}), subtypeStr);

    % Combine all sorting keys
    if meshsz == 1
        meshSize = cellfun(@(x) str2double(x{1}), meshSizeStr);
        sortKeys = [category(:), subtype(:), polyOrder(:), meshSize(:)];
        [~, order] = sortrows(sortKeys, [1 2 3 4]);
    else
        sortKeys = [category(:), subtype(:), polyOrder(:)];
        [~, order] = sortrows(sortKeys, [1 2 3]);
    end

    % Apply sorting
    files_sorted = files(order);
end

