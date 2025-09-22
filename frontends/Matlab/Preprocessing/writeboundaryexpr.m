function writeboundaryexpr(filename, boundaryexpr)

% Define your string array
%boundaryexpr = ["abs(y + 5) < 1e-8", "abs(x - 5) < 1e-8", "abs(y - 5) < 1e-8", "abs(x + 5) < 1e-8"];

% Open a binary file for writing
fid = fopen(filename, 'w');

% Write number of strings
n = numel(boundaryexpr);
fwrite(fid, n, 'int32');

% Write each string: first write its length, then its characters
for i = 1:n
    str = boundaryexpr(i);
    bytes = uint8(char(str));  % Convert string to bytes
    len = length(bytes);
    fwrite(fid, len, 'int32');     % Write string length
    fwrite(fid, bytes, 'uint8');   % Write actual characters
end

fclose(fid);
