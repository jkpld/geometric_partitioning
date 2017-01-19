function groups = triGroups(triList)
% This will recursively find all groups of triangles that are joind by an
% edge.
%
% Input Parameters:
% list - 3 x N array of vertex indices. Each column represents a triangle.
%        (This is similar to what is returned by delaunay -- maybe
%        transposed.)
%
% Output Parameters:
% groups - A cell array with where each element containes a list of
%          triangles that are connected to each other.

% James Kapaldo
% 2016-10-07
% 2016-10-23 : re-wrote triGroup below.

if isempty(triList)
    groups = [];
else
    % Get a triangle group

%     group = triGroup(triList);
    group = triGroup(triList(:,1),triList(:,2:end));
    
    % Remove the group from the list.
    toRemove = sum(group(1,:)' == triList(1,:),1) & sum(group(2,:)' == triList(2,:),1) & sum(group(3,:)' == triList(3,:),1);
    triList = triList(:,~toRemove);

    % Get all of the groups.
    groups = [{group},triGroups(triList)];
    
end

end


function group = triGroup(initialTriangles,triList)
% This will recursively find all triangles that are joined by an edge.

% 2016-20-23 : re-written to have two inputs and no for-loop. no we will
% not get duplicates. old version below.

% Find the indices of all triangles that have to vertices in common with
% the first triangle.
if isempty(triList)
    group = initialTriangles;
    return;
end

pairedTri = false(size(triList,2),1);
for i = 1:size(initialTriangles,2)
    pairedTri = pairedTri | permute(sum(sum(initialTriangles(:,i) == permute(triList,[3,1,2]),1),2)>1,[3,2,1]);
end

if any(pairedTri) 
    % If there are triangles with two vertices in common, then first remove
    % the current triangle from the list and add it to the group
    
    % Add the first triangle to the group and remove it from the list.
    group = [initialTriangles, triGroup(triList(:,pairedTri), triList(:,~pairedTri))];
else
    % If there are no trianges connected to the first triangle than the
    % first triangle is in its own group.
    group = initialTriangles;
end

end


% function group = triGroup(triList)
% % This will recursively find all triangles that are joined by an edge.
% 
% % Find the indices of all triangles that have to vertices in common with
% % the first triangle.
% 
% pairedTri = permute(sum(sum(triList(:,1) == permute(triList(:,2:end),[3,1,2]),1),2)>1,[3,2,1]);
% pairedTriInd = find(pairedTri);
% 
% if ~isempty(pairedTriInd)
%     % If there are triangles with two vertices in common, then first remove
%     % the current triangle from the list and add it to the group
%     
%     % Add the first triangle to the group and remove it from the list.
%     group = triList(:,1);
%     triList(:,1) = [];
% 
%     % Now iterate over each of the connected triangles to find their
%     % connected triangles
%     for i = 1:numel(pairedTriInd)
%         group = [group, triGroup([triList(:,pairedTriInd(i)), triList(:,~pairedTri)])]; %#ok<AGROW>
%     end
% else
%     % If there are no trianges connected to the first triangle than the
%     % first triangle is in its own group.
%     group = triList(:,1);
% end
% 
% end