function vtuwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields)

filename = filename + ".vtu";

% get dimensions
[npoints,nd] = size(cgnodes);
[ncells,nve] = size(cgcells);

if nd==2
    cgnodes = cat(2,cgnodes,zeros(npoints,1));
end
outputCG = zeros(npoints,3);

[~,~,endian] = computer;
if strcmp(endian,'L')
    byte_order = 'LittleEndian';
else
    byte_order = 'BigEndian';
end


formatype = 'appended';

% float and integer byte size
fbytesize = 4; % Float32 format
ibytesize = 4; % Int32   format
obytesize = 8; % Int64   format

nsc = length(scalars);
if nsc>1
    sidx = 1:2:(nsc-1);
    sidy = 2:2:nsc;
else
    sidx = [];
    sidy = [];
end
nvt = length(vectors);
if nvt>1
    vidx = 1:2:(nvt-1);
    vidy = 2:2:nvt;
else
    vidx = [];
    vidy = [];    
end

% Open VTK output file
fid = fopen(filename, 'w');
% VTK DataFile Version
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="%s" header_type="UInt64">\n',byte_order);
fprintf(fid,'  <UnstructuredGrid>\n');
fprintf(fid,'    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',npoints,ncells);

offset = 0;
if (~isempty(sidx)) || (~isempty(vidx))
    fprintf(fid,'      <PointData Scalars="scalars">\n');
end
for ii = 1:length(sidx)
    title = scalars{sidx(ii)};% varargin{sidx(ii)+1};
    fprintf(fid,'        <DataArray type="Float32" Name="%s" Format="%s" offset="%d">\n', title, formatype, offset);
    fprintf(fid,'        </DataArray>\n');
    offset = offset + npoints*fbytesize + obytesize;
end
for ii = 1:length(vidx)
    title = vectors{vidx(ii)};%varargin{vidx(ii)+1};
    fprintf(fid,'        <DataArray type="Float32" Name="%s" NumberOfComponents="3" Format="%s" offset="%d">\n', title, formatype, offset);
    fprintf(fid,'        </DataArray>\n');
    offset = offset + 3*npoints*fbytesize + obytesize;
end
if (~isempty(sidx)) || (~isempty(vidx))
    fprintf(fid,'      </PointData>\n');
end
fprintf(fid,'      <Points>\n');
fprintf(fid,'        <DataArray type="Float32" Name="points" NumberOfComponents="3" Format="%s" offset="%d">\n', formatype, offset);
fprintf(fid,'        </DataArray>\n');
fprintf(fid,'      </Points>\n');
offset = offset + 3*npoints*fbytesize + obytesize;
fprintf(fid,'      <Cells>\n');
fprintf(fid,'        <DataArray type="Int32" Name="connectivity" Format="%s" offset="%d">\n', formatype, offset);
fprintf(fid,'        </DataArray>\n');
offset = offset + ncells*nve*ibytesize + obytesize;
fprintf(fid,'        <DataArray type="Int32" Name="offsets" Format="%s" offset="%d">\n', formatype, offset);
fprintf(fid,'        </DataArray>\n');
offset = offset + ncells*ibytesize + obytesize;
fprintf(fid,'        <DataArray type="UInt8" Name="types" Format="%s" offset="%d">\n', formatype, offset);
fprintf(fid,'        </DataArray>\n');
fprintf(fid,'      </Cells>\n');
fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');
fprintf(fid,'  <AppendedData encoding="raw">\n');
fprintf(fid,'   _');

%  Write all the scalar dataset first
for ii = 1:length(sidx)
    output = fields(:,scalars{sidy(ii)},:);%varargin{ sidx(ii) + 2};
    outputCG(cgelcon(:),1) = output(:);
    fwrite(fid, npoints*fbytesize, 'int64');
    fwrite(fid, outputCG(:,1), 'float32');
end

%  Write all the vector datasets then
for ii = 1:length(vidx)
    tmp = fields(:,vectors{vidy(ii)}(1),:);    
    output = zeros(numel(tmp),3);        
    output(:,1) = tmp(:);
    tmp = fields(:,vectors{vidy(ii)}(2),:);
    output(:,2) = tmp(:);
    if nd==3
        tmp = fields(:,vectors{vidy(ii)}(3),:);
        output(:,3) = tmp(:);        
    end
    % Convert to CG Field
    outputCG(cgelcon(:),1) = output(:,1);
    outputCG(cgelcon(:),2) = output(:,2);
    outputCG(cgelcon(:),3) = output(:,3);
    fwrite(fid, 3*npoints*fbytesize, 'int64');
    fwrite(fid, outputCG', 'float32');
end

output = cgnodes';
fwrite(fid, 3*npoints*fbytesize, 'int64');
fwrite(fid, output, 'float32');

output = cgcells';
fwrite(fid, nve*ncells*ibytesize, 'int64');
fwrite(fid, output, 'int32');

output = nve:nve:nve*ncells;
fwrite(fid, ncells*ibytesize, 'int64');
fwrite(fid, output, 'int32');

output = celltype*ones(1,ncells);
fwrite(fid, ncells, 'int64');
fwrite(fid, output, 'uint8');

fprintf(fid,'  </AppendedData>\n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);

end

