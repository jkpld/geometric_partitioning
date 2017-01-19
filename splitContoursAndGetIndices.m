function [B_cell, n_cell, curvature_cell, cutVertices_cell, cutVertices, ...
          cutVerticesContourNumber, cutCurvatures, B_contourNumber] ...
          = splitContoursAndGetIndices(B, n, curvature, cuts, ncuts)
% splitContourAndGetIndices will split the boundary and the boundary
% normals into it's separate contours which are delimited by rows of NaN's.
% Each seperate contour will be placed in a cell array.
% The cut vertices' boundary index will be determined. A matrix of these
% indices will be output (cutVertices) where the index cooresponds to the
% numbering for each contour in the cell array B_cell. The contour that
% each vertex is in will be returned as cutVerticesContourNumber, which is
% a matrix the same size as cutVertices.
%
% The cutVertices will also be output as a cell array (cutVertices_cell)
% where the vertices from each contour are grouped together and sorted.
% This array will be used later when optimizing unpaired cuts.
% The contour number of each boundary vertex in the input boundary array
% will be output as B_contourNumber.
%
% Input parameters:
%
% B - Mx2 array with the boundary. Seperate contours (outside, holes)
%     should be delimeted by a row of NaN's.
%
% n - Mx2 array with the boundary normals (pointing inward). Seperate
%     counters should be delimted by a row of NaN's.
%
% curvature - Mx1 array with boundary curvature.
%
% cuts - Nx4 array where N is the number of cuts. cuts(i,1:2) and
%        cuts(i,3:4) make up the two vertices of the i'th cut.
%
% ncuts - Nx4 array where N is the number of cuts. ncuts(i,1:2) gives the
%         boundary normal at cuts(i,1:2) and ncuts(i,3:4) gives the
%         boundary normal at cuts(i,3:4).
%
% Output parameters:
%
% B_cell - A cell array containing the boundary for each seperate contour.
%
% n_cell - A cell array containing the boundary normals for each seperate
%          contour.
%
% curvature_cell - A cell array containing the boundary curvatures for each
%                  seperate contour.
%
% cutVertices_cell - A cell array containing the cut vertices in each
%                    contour. The vertices in each cell will be sorted.
%
% cutVertices - A matrix of size N x 2, where N is the number of cuts.
%               cutVertices(i,1) gives the index of cuts(i,1:2). The index
%               returned is for the contour that cutVertices(i,1) is in.
%               This contour number is returned in
%               cutVerticesContourNumber(i,1).
%
% cutVerticesContourNumber - A matrix of size N x 2.
%                            cutVerticesContourNumber(i,1) gives the
%                            contour number (which says what cell of
%                            B_cell, n_cell, and cutVertices_cell the
%                            vertix is appart off) of vertex cuts(i,1:2)
%                            and cutVertices(i,1).
%
% cutCurvatures - Nx4 array with the curvatures of each cut vertex.
%
% B_contourNumber - Mx1 array that gives the contour number for each vertex
%                   in the input array B.
%
% See also CREATEOBJECTPARTITIONS

% James Kapaldo
% 2016-10-06


% First need to find the indices of the cut vertices. Before doing this, we
% need to offset all boundaries and cut vertices by 0.25 in the direction
% of their normal vectors. This is required because holes in the boundary
% may have pixels that overlap. Offseting the boundaries means that no
% pixels should overlap and we can find the proper indices.

offset_B = B + 0.25*n;
offset_cuts = cuts + 0.25*ncuts;

% Find the indices of the cut vertices.
tmpCuts = reshape(offset_cuts',2,size(cuts,1)*2);
[cutVertices,~] =  find( (offset_B(:,1) == tmpCuts(1,:)) & (offset_B(:,2) == tmpCuts(2,:)) );

cutVertices = reshape(cutVertices,2,size(cuts,1))';

% Find the rows of NaN's.
nanLocs = isnan(B(:,1));

% Create the contour numbers.
B_contourNumber = cumsum(nanLocs) + 1;
cutVerticesContourNumber = B_contourNumber(cutVertices);

if size(cutVertices,1) == 1
    cutVerticesContourNumber = cutVerticesContourNumber';
end

% Get vertices' curvature
cutCurvatures = curvature(cutVertices);
if size(cutVertices,1) == 1
    cutCurvatures = cutCurvatures';
end

if ~isempty(nanLocs)
    
    % Initialize cell arrays
    nanLocs = find(nanLocs);
    
    B_cell = cell(numel(nanLocs),1);
    n_cell = cell(numel(nanLocs),1);
    curvature_cell = cell(numel(nanLocs),1);
    cutVertices_cell = cell(numel(nanLocs),1);
    
    % Get the contour bounds
    nanLocs = [0;nanLocs;size(B,1)+1];
    
    for i = 1:numel(nanLocs)-1
        % Save each contour
        B_cell{i} = B(nanLocs(i)+1:nanLocs(i+1)-1,:);
        n_cell{i} = n(nanLocs(i)+1:nanLocs(i+1)-1,:);
        curvature_cell{i} = curvature(nanLocs(i)+1:nanLocs(i+1)-1,:);
        
        % Find the cut vertices in this contour.
        tmpInds = cutVerticesContourNumber == i;
        
        % Offset the cut vertices' indexes for each contour.
        cutVertices(tmpInds) = cutVertices(tmpInds) - nanLocs(i);
        
        % Sort the cut vertices of each contour and save the results.
        cutVertices_cell{i} = sort(cutVertices(tmpInds));
    end
else
    B_cell{1} = B;
    n_cell{1} = n;
    curvature_cell{1} = curvature;
    cutVertices_cell{1} = sort(cutVertices(:));
end

end