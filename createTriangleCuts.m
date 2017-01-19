function [cut,cutN,cutK,cutInd,Info] = createTriangleCuts(B, n, K, triangleGroups, centers, cutGroupAssociation, options)
% CREATETRIANGLECUTS  Create cuts from the from the centers of the
% triangles to the boundary and optimize but the boundry vertex location
% and the triangle center location
% 
% [cut,cutN,cutK,cutInd,Info] = createTriangleCuts(B, n, K, triangleGroups, centers, cutGroupAssociation, options)
%
% Input paramters :
%
% B - Mx2 array with the boundary
%
% n - Mx2 array with the boundary unit normals (pointing inward)
%
% K - Mx1 array with boundary curvature
%
% triangleGroups - Cell array of triangle groups as returned by
%                  constrainedObjectCenterTriangulation
%
% centers - Nx2 array with the object centers
%
% cutGroupAssociation - cell array where each element contains the indices
%                       of the cut vertices associated with the each
%                       triangleGroup. This is optional. If supplied, then
%                       the triangle vertices of group i will be
%                       constrained by the cuts not associated with group i
%                       when the triangle vertices are optimized.
%
% options - element of the class declumpOptions
%
% Output parameters :
%
% cut - Lx4 array where the i'th cut is between the vertices given by
%       cuts(i,1:2) and cuts(i,3:4)
%
% cutN - Lx4 array with the boundary normals of each cut vertex
%
% cutK - Lx2 array with the boundary curvature of each cut vertex
%
% cutInd - Lx2 array with the boundary index of each cut vertex
%
% Info - structure array (1 x numel(triangleGroups)) with the fields
%       triGroup - The triangles in the current group
%       originalCuts - cuts before any optimization
%       vertexOptimizedCuts - cuts after vertex optimizaiton
%       optimizedCuts - cuts after both vertex optimization and triangle
%                       center optimization
%       triangleCenterSearchPoints - the points that were search for the
%                                    optimized triangle center.
%
% See also CONSTRAINEDOBJECTCENTERTRIANGULATION FINDASSOCIATEDCUTS

% James Kapaldo
% 2016-10-19



% Constants :::::::::::::::::::::::::::::::
TRIANGLE_CENTER_SEARCH_FRAC = 0.7;
% :::::::::::::::::::::::::::::::::::::::::


% Preliminaries:
% -------------------------------------------------------------------------

if isempty(triangleGroups)
    cut = {};
    cutN = {};
    cutK = {};
    cutInd = {};
    Info = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
    return;
end


if nargin < 6
    cutGroupAssociation = {};
end

for i = numel(cutGroupAssociation):-1:1
    notAssociated{i} = cat(1,cutGroupAssociation{[1:i-1,i+1:end]});
end

DEBUG = (nargout > 4);

% Get the x-y centers
centers_x = centers(:,1);
centers_y = centers(:,2);

% Get a boundary offset by the inward pointing normal vectors.
offset_B = B + 0.25*n;
%     This offset boundary will be used when finding the index of an
%     intersection with the boundary. We need to offset it because it could
%     be that the boundary of a hole goes through the same pixel twice on
%     opposite sices of the contour. This would give two answers for the
%     the index of intersection at that pixel. Offseting the boundary
%     ensures that there are no repeated x-y values, and so the index of
%     intersection will be unique.


% Initialize arrays to store any triangle cuts
cut = cell(1,numel(triangleGroups));
cutN = cell(1,numel(triangleGroups));
cutK = cell(1,numel(triangleGroups));
cutInd = cell(1,numel(triangleGroups));


if DEBUG
    Info(numel(triangleGroups)) = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
end

% -------------------------------------------------------------------------
% End preliminaries

for groupNum = 1:numel(triangleGroups)
%     groupNum
    % Extract the current group of triangles. Remember that if this group
    % contains more than one triangle (size(triGroup,2)>1) then the
    % triangles share an edge.
    triGroup = triangleGroups{groupNum}; % 3xQ
    
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
    
    % Find the shared edges of the triangles.
    [~,~,ic1] = unique(triEdges,'rows');
    doubleEdges = accumarray(ic1,1)>1;
    sharedEdges = sum((ic1 == find(doubleEdges')) .* (1:sum(doubleEdges)), 2); % array of size 3*Q x 1, each pair will have a different number.

    % Get the centers of the triangles
    triCenters = [ mean(centers_x(triGroup),1)', mean(centers_y(triGroup),1)' ]; % Qx2
    
    % Get the mid point of the triangle edges
    edgeMidPoints = [ mean(centers_x(triEdges),2), mean(centers_y(triEdges),2) ]; % 3*Qx2
    
    % Get the triangle edge normals. These edge normals will not all be
    % pointing outward, or all point inward, but will be mixed.
    edgeNormals = centers(triEdges(:,1),:) - centers(triEdges(:,2),:);
    edgeNormals = [-edgeNormals(:,2),edgeNormals(:,1)];
    edgeNormals = edgeNormals ./ sqrt(sum(edgeNormals.^2,2)); % 3*Qx2
    
    % Get the unit vectors from the triangle centers to the edgeMidPoints
    uv = edgeMidPoints - kron(triCenters,ones(3,1));
    uv = uv ./ sqrt(sum(uv.^2,2)); % 3*Qx2
    
    % Correct the edge normals to all point outward
    edgeNormals = sign(sum(uv.*edgeNormals,2)) .* edgeNormals;
    
    if DEBUG
        Info(groupNum).triGroup = triGroup;
    end
    
    
    % Create new potential cuts
    % ....................................................................
    %
    % Iterate over each triangle in the group.
    % * Create new potential cuts from the center of the triangle through
    %   the mid point of any un-shared triEdges to the nearest boundary.
    % * Optimize the new potential cuts through un-shared triEdges.
    % * Create a temporary cut from the center of the triangle to the mid
    %   point of shared triEdges.
    % * Optimize the center of the triangle - do not consider the dot
    %   product between the temporary cut and the shared edge normal.

    numCuts = numel(sharedEdges) - sum(logical(sharedEdges))/2;
    newPotentialCuts = zeros(numCuts,4);
    newPotentialCutNs = nan(numCuts,4);
    newPotentialCutKs = nan(numCuts,2);
    newPotentialCutInds = nan(numCuts,1);
    
    if DEBUG
        originalCuts = zeros(sum(~sharedEdges),4);
        vertexOptimizedCuts = zeros(sum(~sharedEdges),4);
    end
    
    counter = 1;
%     triGroup
    
%     figure(55)
%         clf(55)
    
    for triNum = 1:size(triGroup,2)
%         triNum
        cntInds = (1:3) + (triNum-1)*3;
        isShared = sharedEdges(cntInds)>0;
        
        % Part 1: 
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Find the closest boundary vertex to the triangle edge mid points.
        % Only consider boundary vertices that have a positive dot product
        % with the outward pointing triange edge normal.

        tmpCounter = 1;
        pCutInds = zeros(sum(~isShared),1);
        for j = cntInds(~isShared)
            tmpB = offset_B - edgeMidPoints(j,:);
            possibleInds = find(sum(edgeNormals(j,:) .* tmpB,2) > 0);
            [~,closestInd] = min(sum(tmpB(possibleInds,:).^2,2));
            pCutInds(tmpCounter) = possibleInds(closestInd);
            tmpCounter = tmpCounter + 1;
        end
        
%         figure(111)
%         clf(111)
%         line(B(:,1),B(:,2),'color','w')
%         line(triCenters(triNum,1),triCenters(triNum,2),'marker','x','markerSize',15,'color','r')
%         line(centers_x,centers_y,'linestyle','none','marker','*','color','g')
%         hold on
%         quiver(edgeMidPoints(cntInds(~isShared),1),edgeMidPoints(cntInds(~isShared),2),edgeNormals(cntInds(~isShared),1),edgeNormals(cntInds(~isShared),2),'color','r')
%         daspect([1 1 1])
%         drawnow
%         goDark(gcf)
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % End Part 1: 
        
        numPCuts = numel(pCutInds);
        
        if DEBUG
            originalCuts(counter:counter+numPCuts-1,:) = [B(pCutInds,:), triCenters(triNum,:) .* ones(numPCuts,1)];
        end
        
        originalpCutVrts = B(pCutInds,:);
        % Part 2:
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Optimize the new potential cuts through un-shared triEdges.
        
        
        walls = pCutInds;
        if ~isempty(notAssociated)
            walls = sort([walls(:);notAssociated{groupNum}]);
        end
        
        [pCutVrts,pCutVrtNs,pCutVrtKs,pCutInds] = optimizeVertex(pCutInds,triCenters(triNum,:),walls,B,n,K,options);
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % End Part 2
        
        if DEBUG
            vertexOptimizedCuts(counter:counter+numPCuts-1,:) = [pCutVrts, triCenters(triNum,:) .* ones(numPCuts,1)];
        end
        
        
        % Part 3:
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Create a temporary cut from the center of the triangle to the mid
        % point of shared triEdges.
        tmpCutVrts = edgeMidPoints(cntInds(isShared),:);
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % End Part 3
        
        
        
        % Part 4:
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Optimize the center of the triangle - do not consider the dot
        % product between the temporary cut and the shared edge normal.
        
        % Optimize the center of the triangle over the points inside of
        % a second triangle whose vertices are a fraction along the
        % exsisting cuts.
        triangleSearchBounds = (1-TRIANGLE_CENTER_SEARCH_FRAC)*triCenters(triNum,:) + TRIANGLE_CENTER_SEARCH_FRAC * [pCutVrts;tmpCutVrts];
        minSearchBounds = round(min(triangleSearchBounds) - 1);
        maxSearchBounds = round(max(triangleSearchBounds) + 1);
        
        % Create the x and y values to search over such that the
        % triangle center is always one of the search values.
        x = [triCenters(triNum,1):-2:minSearchBounds(1), triCenters(triNum,1)+2:2:maxSearchBounds(1)];
        y = [triCenters(triNum,2):-2:minSearchBounds(2), triCenters(triNum,2)+2:2:maxSearchBounds(2)];
        X = x' .* ones(1,numel(y));
        Y = y  .* ones(numel(x),1);
%         X
%         Y
        search_r = [X(:),Y(:)];
%         search_r
%         triangleSearchBounds(:,1)
%         triangleSearchBounds(:,2)
%         minSearchBounds
%         maxSearchBounds
        
        
%         line(B(:,1),B(:,2),'color','w')
%         line(triCenters(triNum,1),triCenters(triNum,2),'marker','x','markerSize',15,'color','r')
%         line(centers_x,centers_y,'linestyle','none','marker','*','color','g')
%         line(originalpCutVrts(:,1),originalpCutVrts(:,2),'linestyle','none','marker','*','color','m')
%         line(pCutVrts(:,1),pCutVrts(:,2),'linestyle','none','marker','*','color','c')
% %         line(triangleSearchBounds(:,1),triangleSearchBounds(:,2),'color','r','marker','o')
%         drawnow;
%         goDark(gcf);
%         
        % Remove any search point that is outside of the search
        % triangle.
        in = inpolygon(search_r(:,1),search_r(:,2),triangleSearchBounds(:,1),triangleSearchBounds(:,2));
        search_r(~in,:) = [];
        
%         line(search_r(:,1),search_r(:,2),'marker','s','linestyle','none','markersize',5,'color','y')
        
        if ~isempty(search_r)
            
            r1 = search_r - permute(pCutVrts,[3,2,1]); % Lx2xK
            d1 = sqrt(sum(r1.^2,2)); % Lx1xK
            r1 = r1 ./ d1; % Lx2xK

            if any(isShared)
                r2 = search_r - permute(tmpCutVrts,[3,2,1]); % Lx2xR
                d2 = sqrt(sum(r2.^2,2)); % Lx1xR
                d1 = cat(3,d1,d2); % Lx1x(K+R);
            end

            d1 = mean(d1,3);

            a1 = mean(sum(r1 .* permute(pCutVrtNs,[3,2,1]),2),3);

            valueToOptimize = a1 ./ d1;

            [~,bestTriCentInd] = max(valueToOptimize);

            triCenters(triNum,:) = search_r(bestTriCentInd,:);
        end
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % End Part 4
        
%         line(triCenters(triNum,1),triCenters(triNum,2),'marker','o','markerSize',15,'color','g')
%         
%         for j = 1:numPCuts
%             line([pCutVrts(j,1), triCenters(triNum,1)], [pCutVrts(j,2), triCenters(triNum,2)],'Color','c','LineWidth',2)
%         end
%         daspect([1 1 1])
        
        if DEBUG
            Info(groupNum).triangleCenterSearchPoints = [Info(groupNum).triangleCenterSearchPoints; search_r];
        end
        
        % Save the results for later.
        triCenters(triNum,:) = round(triCenters(triNum,:));
        newPotentialCuts(counter:counter+numPCuts-1,:) = [pCutVrts, triCenters(triNum,:) .* ones(numPCuts,1)];
        newPotentialCutNs(counter:counter+numPCuts-1,1:2) = pCutVrtNs;
        newPotentialCutKs(counter:counter+numPCuts-1,1) = pCutVrtKs;
        newPotentialCutInds(counter:counter+numPCuts-1,:) = pCutInds;
        
        counter = counter + numPCuts;
    end % for each triangle
    
    if DEBUG
        Info(groupNum).originalCuts = originalCuts;
        Info(groupNum).vertexOptimizedCuts = vertexOptimizedCuts;
    end
       
    % ....................................................................
    % * Create new potenial cuts between each of the triangle centers that
    %   share an edge

    for sharedEdgeNum = 1:max(sharedEdges)
        if counter > size(newPotentialCuts,1)
            warning('size allocation wrong')
        end
        pairedCenters = ceil(find(sharedEdges==sharedEdgeNum)/3);
        newPotentialCuts(counter,:) = [triCenters(pairedCenters(1),:), triCenters(pairedCenters(2),:)];
        counter = counter + 1;
        % Create a cut between the two pairedCenters.
    end
    % ....................................................................
    
    cut{groupNum} = newPotentialCuts;
    cutN{groupNum} = newPotentialCutNs;
    cutK{groupNum} = newPotentialCutKs;
    cutInd{groupNum} = newPotentialCutInds;
    
    
    if DEBUG
        Info(groupNum).optimizedCuts = newPotentialCuts;
    end    
end


end





function [cut,cutN,cutK,optimizedVertexIndices] = optimizeVertex(vertexIndices,center,wallIndices,B,n,K,options)
% OPTIMIZEVERTEX  Search in the local region of a vertex for a new vertex
% that minimizes the distance to center, maximizes the dot product between
% the vertex normal and the vertex-center direction, and maximizes the
% vertex curvature.
%
% Input parameters:
%
% vertexIndices - Nx1 array with the indices of the vertices to optimize.
%
% B - Mx2 array with the boundary
%
% n - Mx2 array with the boundary unit normals (pointing inward)
%
% K - Mx1 array with boundary curvature
%
% wallIndices - Lx1 array of wall indices. vertexIndices must be a subset
%               of this array.
%
% options - structure with one required field
%           searchRadius - The radius of the search range over which the
%           unpaired cut will be optimized.
%
% Output parameters:
%
% vertex - Nx2 array where vertex(i,:) is the optimized vertex for
%          vertexIndices(i)
%
% vertexN - Nx2 array where vertexN(i,:) is the unit normal at the
%           optimized vertex vertex(i,:)
%
% vertexK - Nx1 array where vertexK is the curvature at the optimized
%           vertex(i,:)
%
% optimizedVertexIndices - Nx1 array with the optimized vertex indices.
%
% A vertex is optimized by searching for a new vertices within searchRadius
% of the current vertices such that the vertex-center distance is
% minimized, the dot product of the cut' vertex normals with the
% vertex-center direction is maximized, and the vertex curvature is
% maximized. The search region is not allowed to cross any other cut vertex
% (which would cause cuts to potentially cross each other).
%
% See also GETSEARCHINDS CREATEOBJECTPARTITIONS

% James Kapaldo
% 2016-10-16
% 2016-10-19 : removed positive curvature limit.

% Getting the search range around each unpairedCut vertex requires some
% care. We cannot simply use B(i-searchRange:i+searchRange), but the search
% range must 1) wrap around each contour, and 2) only contain indices of
% the same contour as the upairedCut vertex. Furthermore, the searchRange
% should not cross any other cut vertex.

% Get search range.
SEARCH_RADIUS = options.Search_Radius;



% Put limits on the curvature;
WIGNER_SEITZ_RADIUS = options.Wigner_Seitz_Radius;
maxK = 1/WIGNER_SEITZ_RADIUS;
% K(K>maxK) = maxK;
K(K<-maxK) = -maxK;
% K(K<0) = K(K<0)*5;


% Get the contours of the boundary
nanLocs = isnan(B(:,1));
BcontourNumber = cumsum(nanLocs) + 1;
Binds = (1:size(B,1))';

vertexContourNumbers = BcontourNumber(vertexIndices);
wallContourNumbers = BcontourNumber(wallIndices);


walls = accumarray(wallContourNumbers(:), wallIndices(:),[sum(nanLocs)+1,1],@(x) {sort(x)},[]);
contourInds = accumarray(BcontourNumber(~nanLocs),Binds(~nanLocs),[sum(nanLocs)+1,1],@(x) {x},[]);

optimizedVertexIndices = zeros(size(vertexIndices,1),1);

for i = 1:size(vertexIndices,1)

    % Get the search indices for the two vertices of the current cut
    c1 = vertexContourNumbers(i);
    [searchInds1, ind1] = getSearchInds(vertexIndices(i),contourInds{c1},walls{c1},SEARCH_RADIUS);

    % Get the locations and normals over each sarch range.
    r1 = B(searchInds1,:); % N1 x 2
    n1 = n(searchInds1,:);
    k1 = K(searchInds1,:);

    % Get the vector between each possible search location
    r = r1 - center; % N1 x 2
    d = sqrt(sum(r.^2,2)); % N1 x 1
    r = r ./ d; % N1 x 2
    
    % Get the dot products between the vertex normals and the cut
    % directions
    a1 = sum( -r .* n1, 2); % N1 x 1
    
    % Optimize
    valueToOptimize = (a1+k1) ./ d; % N1 x 1
    [~,ind] = max(valueToOptimize);
    
    % Update the wall positions.
    walls{c1}(ind1) = searchInds1(ind);

    % Save the new indices
    optimizedVertexIndices(i) = searchInds1(ind);
    
end

cut = B(optimizedVertexIndices,:);
cutN = n(optimizedVertexIndices,:);
cutK = K(optimizedVertexIndices);

end