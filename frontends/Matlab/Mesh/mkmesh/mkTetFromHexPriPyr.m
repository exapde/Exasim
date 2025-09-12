function [ t ] = mkTetFromHexPriPyr(t)
%MKTETFROMHEXPRIPYR - Create Tetrahedral mesh from Hexa/Prism/Pyramid mesh.
% This approach modifies the table of connectivities to split non tets
% (i.e. hexahedra, prisms or pyramids) into tets only. This approach is
% based only on the table of connectivities and not the nodes coordinates.
% We follow the document "How to subdivide Pyramids, Prisms ans Hexahedra
% into Tetrahedra", from J. Dompierre, P. Labbe, M.-G. Vallet and R. 
% Camarero, Rapport CERCA R99-78, 8th International Meshing Roundtable,
% 1999.
%
% ATTENTTION : Currently, only the Hexa can be broken with the current
% version. Prism and Pyramids have to be added, althought the theory is
% ready. NON-CONFORMING MESH MAY NOT WORK.
%
% SYNTAX:  [t] = mkTetFromHexPriPyr(p,t)
%
% INPUTS:
%    t - Elements connectivity table
%
% OUTPUTS:
%    t - Elements connectivity table, including tets only.
%
%
% OTHER M-FILES REQUIRED: none
% SUBFUNCTIONS: none
% MAT-FILES REQUIRED: none
%
% SEE ALSO: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Author(s): Terrana Sebastien
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email address: terrana@mit.edu 
% Website: http://aeroastro.mit.edu/
% January 2018

% Gets dimensions
[ne, nv] = size(t);


if nv == 5 % Pyramids
    % Initialize Tets (2 tets per pyramid)
    t2 = zeros(2*ne,4);
    % Minimum diagonal for each pyramid
    mindiag = min(t(:,[1 3]),[],2) < min(t(:,[2 4]),[],2);
    % Local Indices of new Tets
    IndTets = [2 3 4 5 2 4 1 5; ...
               1 2 3 5 1 3 4 5];
    % Create the New Tets 2 by 2
    for i=1:ne
        t2(2*i-1,:)= t(i,IndTets(mindiag(i)+1,1:4));
        t2(2*i,:)  = t(i,IndTets(mindiag(i)+1,5:8));
    end
    % Create output Tets table
    t = t2;
    
elseif nv == 6 % Prisms
    % Initialize Tets (3 tets per prism)
    t2 = zeros(3*ne,4);
    % Gets the smallest vertex for each prism
    [~,indminv] = min(t,[],2);
    % Indirections of vertices according to the smallest vertex
    IndVert = [1 2 3 4 5 6; 2 3 1 5 6 4; 3 1 2 6 4 5; ...
               4 6 5 1 3 2; 5 4 6 2 1 3; 6 5 4 3 2 1];
    % Create the new tets 3 by 3
    for i=1:ne
        % Get the indirection for the rotated prism
        indPrism = t(i,IndVert(indminv(i),:));
        % Split the rorated prism into 3 tets, there are 2 ways to do that,
        % according to how the third quadrilateral face is split.
        if min(indPrism([2 6]))<min(indPrism([3 5]))
            t2(3*i-2,:)= indPrism([1 2 3 6]);
            t2(3*i-1,:)= indPrism([1 2 6 5]);
            t2(3*i,:)  = indPrism([1 5 6 4]);
        else
            t2(3*i-2,:)= indPrism([1 2 3 5]);
            t2(3*i-1,:)= indPrism([1 5 3 6]);
            t2(3*i,:)  = indPrism([1 5 6 4]);
        end
    end
    % Create output Tets table
    t = t2;
        
elseif nv == 8 % Hexahedra
    % Initialize Tets (6 tets per hex)
    t2 = zeros(6*ne,4);
    % Gets the smallest vertex for each prism
    [~,indminv] = min(t,[],2);
    % Indirections of vertices according to the smallest vertex
    IndVert = [1 2 3 4 5 6 7 8; 2 1 5 6 3 4 8 7; 3 2 6 7 4 1 5 8; 4 1 2 3 8 5 6 7; ...
               5 1 4 8 6 2 3 7; 6 2 1 5 7 3 4 8; 7 3 2 6 8 4 1 5; 8 4 3 7 5 1 2 6];
    
    ntet = 0;
    for i=1:ne
        % Get the indirection for the rotated hex
        indHex = t(i,IndVert(indminv(i),:));
        % Gets which diagonals pass through the vertex indHex(7)
        diagBits = [min(indHex([2,7]))<min(indHex([3,6])),...
                    min(indHex([4,7]))<min(indHex([3,8])),...
                    min(indHex([5,7]))<min(indHex([6,8]))];
        % Number of diagonals passing through indHex(7)
        ndiag = sum(diagBits);
        % Second hex rotation around the axis (indHex(1),indHex(7))
        if (all(diagBits-[0 0 1]) || all(diagBits-[1 1 0]))
            % A 120 deg rotation is performed
            indHex = indHex([1 5 6 2 4 8 7 3]);
        elseif (all(diagBits-[0 1 0]) || all(diagBits-[1 0 1]))
            % A 240 deg rotation is performed
            indHex = indHex([1 4 8 5 2 3 7 6]);
        end
        % Split the current hex into tets
        ttets = divideHex2Tets(indHex,ndiag);
        nnewt = size(ttets,1);
        % Add the new tets to the table t2
        t2(ntet+1:ntet+nnewt,:) = ttets;
        % Update total number of tets
        ntet = ntet + nnewt;
    end
    % Create final output table of tets
    t = t2(1:ntet,:);

else
    error('The connectivity table has an unexpected shape.');
end

end


function [ttet] = divideHex2Tets(hex,ndiag)

    switch (ndiag)
        case 0  % Only 5 tets are built in this case
            ttet(1,:) = [hex(1) hex(2) hex(3) hex(6)];
            ttet(2,:) = [hex(1) hex(3) hex(8) hex(6)];
            ttet(3,:) = [hex(1) hex(3) hex(4) hex(8)];
            ttet(4,:) = [hex(1) hex(6) hex(8) hex(5)];
            ttet(5,:) = [hex(3) hex(8) hex(6) hex(7)];
        case 1 % Remark : there is another correct way to spit into tets
            ttet(1,:) = [hex(1) hex(6) hex(8) hex(5)];
            ttet(2,:) = [hex(1) hex(2) hex(8) hex(6)];
            ttet(3,:) = [hex(2) hex(7) hex(8) hex(6)];
            ttet(4,:) = [hex(1) hex(8) hex(3) hex(4)];
            ttet(5,:) = [hex(1) hex(8) hex(2) hex(3)];
            ttet(6,:) = [hex(2) hex(8) hex(7) hex(3)];
        case 2 % Remark : there is another correct way to spit into tets
            ttet(1,:) = [hex(1) hex(5) hex(6) hex(7)];
            ttet(2,:) = [hex(1) hex(4) hex(8) hex(7)];
            ttet(3,:) = [hex(1) hex(8) hex(5) hex(7)];
            ttet(4,:) = [hex(1) hex(2) hex(3) hex(6)];
            ttet(5,:) = [hex(1) hex(4) hex(7) hex(3)];
            ttet(6,:) = [hex(1) hex(7) hex(6) hex(3)];
        case 3 % Remark : there is another correct way to spit into tets
            ttet(1,:) = [hex(1) hex(3) hex(4) hex(7)];
            ttet(2,:) = [hex(1) hex(4) hex(8) hex(7)];
            ttet(3,:) = [hex(1) hex(8) hex(5) hex(7)];
            ttet(4,:) = [hex(1) hex(6) hex(7) hex(5)];
            ttet(5,:) = [hex(2) hex(6) hex(7) hex(1)];
            ttet(6,:) = [hex(2) hex(7) hex(3) hex(1)];
        otherwise
            error('This Diagonal Configuration cannot exist');
    end
end


