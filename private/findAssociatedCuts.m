function associatedCuts = findAssociatedCuts(triangleGroups, cuts, centers, vertexAdjacency, cutCenters)
% FINDASSOCIATEDCUTS  Determine the cuts associated with each triangle
% group and return a cell array of the associated cut indices.
%
% associatedCuts = findAssociatedCuts(triangleGroups, cuts, centers)
%
% triangleGroups - Cell array of triangle groups as returned by
% constrainedObjectCenterTriangulation
%
% cuts - Mx4 array of cuts. cuts(i,1:2) and cuts(i,3:4) give the two
% vertices of the i'th cut.
%
% centers - Nx2 array of object centers
%
% associatedCuts - Cell array of the same size as triangleGroups.
% associatedCuts{i} containes of the indices of cuts associated with the
% triangle group in triangleGroups{i}.
%
% See also CONSTRAINEDOBJECTCENTERTRIANGULATION CREATEOBJECTPARTITIONS

% JamesKapaldo
% 2016-10-19

associatedCuts = cell(size(triangleGroups));
if isempty(triangleGroups)
    return;
end

centers_x = centers(:,1);
centers_y = centers(:,2);

% Get the cut line segments
cutLines = reshape([cuts, nan(size(cuts,1),2)]',2,size(cuts,1)*3)';

for groupNum = 1:numel(triangleGroups)
    
    % Extract the current group of triangles. Remember that if this group
    % contains more than one triangle (size(triGroup,2)>1) then the
    % triangles share an edge.
    triGroup = triangleGroups{groupNum}; % 3xQ
    uniqueTriGroup = unique(triGroup(:));
    
    % Get the edges for each triangle. 
    % * Form the edges (the 2 connected centers) along the 3rd dimension.
    %    ->   cat(3, triGroup, circshift(triGroup,1,1)) 
    %         size : 3xQx2, where N is the number of triangles in the
    %         group
    % * Permute the array to have size 2x3xN
    % * Reshape to have size 2x3*Q
    % * sort so that the smalled edge center is first
    % * Transpose to give final array with size 3*Q x 2, where each row
    %   gives the two centers that make up an edge of the triangle
    triEdges = sort( reshape( permute( cat(3,triGroup,circshift(triGroup,1,1)), [3,1,2]), 2, 3*size(triGroup,2) ))'; % The edges of the i'th triangle are triEdges((1:3) + (i-1)*3,:)

    % Remove duplicate edges
    triEdges = unique(triEdges,'rows');

    % Find the cuts that intersect with the triangle edges
    edges = nan(3*size(triEdges,1),2);
    edges(:,1) = reshape([centers_x(triEdges),nan(size(triEdges,1),1)]',3*size(triEdges,1),1);
    edges(:,2) = reshape([centers_y(triEdges),nan(size(triEdges,1),1)]',3*size(triEdges,1),1);

    [~,~,cutIdx] = intersections(cutLines(:,1),cutLines(:,2),edges(:,1),edges(:,2),1);
    cutIdx = unique(ceil(cutIdx/3));
    
    vrtIdx = reshape(2*cutIdx + [-1,0], 2*numel(cutIdx), 1);
    associated = [cutIdx; ...
                 find(all(reshape(any(vertexAdjacency(:,vrtIdx),2)',2,size(cuts,1))',2))]; % Find other cuts that have both vertices associated with one of the cuts intersected by the triangle edges.
    
    % Also grab in any paired edge that has both centers in the current
    % group. 
    currentCuts = sum(sum(cutCenters == permute(uniqueTriGroup,[3,2,1]),3),2);
	unpairedCutsIdx = find(currentCuts==1 & isnan(cutCenters(:,2)));
    pairedCutsIdx = find(currentCuts==2);
    pairedCuts = cutCenters(pairedCutsIdx,:);
    unpairedCuts = cutCenters(unpairedCutsIdx,1);
    associatedCutCnts = cutCenters(associated,:);
    
    unpairedCutsIdx( ismember(unpairedCuts, [pairedCuts(:);associatedCutCnts(:)] )) = [];

    cutIdx = unique([associated;pairedCutsIdx; unpairedCutsIdx]);
    vrtIdx = reshape(2*cutIdx + [-1,0], 2*numel(cutIdx), 1);
    associated = [cutIdx; ...
                 find(all(reshape(any(vertexAdjacency(:,vrtIdx),2)',2,size(cuts,1))',2))]; % Find other cuts that have both vertices associated with one of the cuts intersected by the triangle edges.
             
    associatedCuts{groupNum} = associated;
end % for



end
% findAssociatedCuts : changeLog
% 2016-10-20 : added cutAdjacency input. if two associated edges are
%              adjacent to both ends of another cut, then that cut is also
%              associated.
% 2016-10-21 : fixed cutadjacency associated to work with sum>2 instead of
%              all(). added in cutCenters to input. if there is only one
%              unpaired cut with a given center, than that cut will be
%              associated with the triangle group (if any) with that
%              center.
% 2016-10-25 : fixed correction of 21st. If there is only one unpaired cut
%              and there is not a paired cut for a center that is a part of
%              a triangle, then associate the unpaired cuts.
% 2016-10-27 : changed the function from working with cut adjacnecies to
%              working with vertex adjacencies. A cut is not adjacent if
%              both of its vertices are also adjacent to any cut that is
%              intersected by the triangle edges or that has both of its
%              center associates in the triangle group or if it is the only
%              edge (unpaired) associated with a center in the group.