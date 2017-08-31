function triangleGroups = constrainedObjectCenterTriangulation(B,centers,options)
% CONSTRAINEDOBJECTCENTERTRIANGULATION  Compute constrained triangulation
% on object centers and return groups of edge-connected triangles in cell
% array.
%
% triangleGroups = constrainedObjectCenterTriangulation(B,centers)
%
% B - Mx2 array with boundary contours. Object and hole contours should be
%     delimited by a row of NaNs.
%
% N - Mx2 array with boundary contours unit normals (pointing inward).
%
% centers - Nx2 array of object centers
%
% options - element of class declumpOptions
%
% triangleGroups - cell array containing arrays of size 3xL where each
% column gives the indices of three centers that form a triangle. All
% triangles returned in a single cell element share edges.
%
% The delaunay triangulation of the centers is computed. If an edge of the
% triangulation intersects with the object boundary, then any triangle with
% that edge is removed. Any triangle as an interior angle larger then set
% threshold (MAX_ANGLE = 135 deg) is removed. Triangle groups are computed
% and returned.
%
% See also COMPUTEOBJECTCENTERS

% James Kapaldo
% 2016-10-19


% Constants :::::::::::::::::::::::::::::::
MAX_ANGLE = options.Max_Interior_Angle * (pi/180);
MIN_ANGLE = options.Min_Interior_Angle * (pi/180);
% :::::::::::::::::::::::::::::::::::::::::

% If there are not three centers then there is no triangle - return.
if size(centers,1) < 3
    triangleGroups = {};
    return;
end

% Get list of triangles that do not intersect the boundary.
% B = B + 3*N;
triangleList = nonBoundaryIntersectingTriangles(B,centers,options.Use_GPU); % Lx3 

% If there are no triangle that do not intersect the boundary, then return.
if isempty(triangleList)
    triangleGroups = {};
    return;
end

% Find groups of triangles. A group of triangles are triangles that share
% an edge.

triangleGroups = triGroups(triangleList'); % cell array with elements of size of 3xLi where sum(Li) = L


% Search through each triangle group for triangles with internal angles
% larger than MAX_ANGLE and remove them. There is one exeption to this: 
%   If three triangles each shair two edges and all have a common vertex,
%   then these three triangles form one larger triangle. Each of these
%   three triangles will have an angle around 120 (if the triangle is
%   equilateral and the middle vertex is at the center of the three outer
%   vertices). Keep these three triangles.

groupsToRemove = false(1,numel(triangleGroups));

for group = 1:numel(triangleGroups)
    

    triangleList = triangleGroups{group};

    toRemove = false(1,size(triangleList,2));
    exempt = false(1,size(triangleList,2));
    
%     % Find exempt triangles
%     if size(triangleList,2) > 2
%         % Look for a large triangle made from three smaller triangles.
%         %  To look for : 1 common vertex for all three triangles. each of
%         %  the other two vertics are apart of another triangle. thus, the
%         %  histogram of the vertex indices of these three triangles should
%         %  be [2 2 2 3]
%         % Start finding subsets with a common vertex
%         
%         vrtList = unique(triangleList(:));
%         A = sum(vrtList == permute(triangleList,[3,2,1]),3);
%         threeCommon = find(sum(A,2)>=3);
%         if ~isempty(threeCommon)
%             for vrt = 1:numel(threeCommon)
%                 currentTris = A(:,logical(A(threeCommon(vrt),:)));
%                 possibleSecondaryVrts = find(sum(currentTris,2)==2);
%                 if numel(possibleSecondaryVrts) == 3
%                     exempt = exempt | sum(A([threeCommon(vrt);possibleSecondaryVrts],:),1)==3;
%                 end
%             end
%         end
% 
%     end
    
    % Find triangles with an internal angle larger than MAX_ANGLE
    for i = 1:size(triangleList,2)
        th = interiorAngles(centers(triangleList(:,i),:));
        if any(th>MAX_ANGLE) || any(th<MIN_ANGLE)
            toRemove(i) = true;
        end
    end
    
    % Remove any triangles that have an internal angle larger than
    % MAX_ANGLE and are not exempt
    triangleList(:,toRemove & ~exempt) = [];
    
    if isempty(triangleList)
        groupsToRemove(group) = true;
    end
    
    triangleGroups{group} = triangleList;
    
end

triangleGroups(groupsToRemove) = [];

end
% constrainedObjectCenterTriangulation : changeLog
% 2016-10-21 : added in the boundary normal to inputs. now we will make
%              sure the triagulation edges to not intersect with the
%              boundary moved in by 3 pixels.
% 2016-10-23 : removed the fix from 2016-10-21, max angle dropped to 110
%              from 135. added in code to allow three triangles that form
%              one larger triangle even if they have an internal angle
%              larger than MAX_ANGLE, but it did not help significantly, so
%              it is commented out.
% 2016-10-29 : added in options input


function triangleList = nonBoundaryIntersectingTriangles(B,centers,useGPU)
% NONBOUNDARYINTERSECTINGTRIANGLES  Compute a list of traingles with
% vertices given by CENTERS whose edges do not overlap with the boundary B.

% James Kapaldo

% Get the daulaunay triangulation
try
    triangleList = delaunay(centers);
catch ME
    if strcmp(ME.identifier, 'MATLAB:delaunay:EmptyDelaunay2DErrId')
        triangleList = [];
        return;
    else
        rethrow(ME)
    end
end

% Get the indices for each edge of the triangulation
triEdges = [triangleList(:,[1,2]); triangleList(:,[2,3]); triangleList(:,[3,1])];
triEdges = unique(sort(triEdges, 2, 'ascend'),'rows');

% From the line segments of each edge as a NaN delimited list.
xVerts = [centers(triEdges(:,1),1) centers(triEdges(:,2),1) nan(size(triEdges,1),1)]';
yVerts = [centers(triEdges(:,1),2) centers(triEdges(:,2),2) nan(size(triEdges,1),1)]';

% Get the intersections of the triangulation edges with the boundary.
[~,~,intrsctns] = intersections(xVerts(:),yVerts(:),B(:,1),B(:,2),0,0,useGPU);

% Get the indices of each edge that intersects with the boundary.
badEdges = triEdges(unique(ceil(intrsctns/3)),:);

% Remove all triangles with an edge that intersects the boundary.
for i = 1:size(badEdges,1)
    triangleList(any(triangleList == badEdges(i,1),2) & any(triangleList == badEdges(i,2),2),:) = [];
end

end
% nonBoundaryIntersectingTriangles : changeLog
% 2016-10-30: added fix for delaunay triangulation when points may be
%             collinear.