function [triCuts, cutsToRemove, Info] = optimizeTriangles( B, n, K, centers, cutCenters, cutIndices, S, options)




% 2016-10-17
% 2016-10-18 : undocumented changes
% 2016-10-19 : Fixed bug in detecting intersections between new triangle
%              cuts and previous cuts. 


% Constants :::::::::::::::::::::::::::::::
MAX_ANGLE = 135 * pi/180;
MAX_CUT_LENGTH = sqrt(sum(range(B).^2));
TRIANGLE_CENTER_SEARCH_FRAC = 0.7;
% :::::::::::::::::::::::::::::::::::::::::

% If there are not three centers then there is no triangle - return.
if size(centers,1) < 3
    triCuts = [];
    cutsToRemove = false(size(cutIndices,1));
    Info = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
    return;
end

% Get list of triangles that do not intersect the boundary

triangleList = nonBoundaryIntersectingTriangles(B,centers); % Lx3 

% If there are no triangle that do not intersect the boundary, then return.
if isempty(triangleList)
    triCuts = [];
    cutsToRemove = false(size(cutIndices,1));
    Info = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
    return;
end

% Remove any triangles that have an internal angle larger than MAX_ANGLE
toRemove = false(size(triangleList,1),1);

for i = 1:size(triangleList,1)
    th = interiorAngles(centers(triangleList(i,:),:));
    if any(th>MAX_ANGLE)
        toRemove(i) = true;
    end
end
triangleList(toRemove,:) = [];

% If there are no triangles left, then return.
if isempty(triangleList)
    triCuts = [];
    cutsToRemove = false(size(cutIndices,1));
    Info = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
    return;
end

% Preliminaries:
% -------------------------------------------------------------------------
% If we reach this point in the code, then computation will be required, so
% we can start defining constants and calculating quantities that will be
% needed later.

DEBUG = (nargout > 2);

% Get cuts, cut normals, and cut curvatures. -- Note we could just pass
% these to the function, but I think fewer function inputs is nice and
% these are very easy to get from the cutIndices.
cuts = [B(cutIndices(:,1),:),B(cutIndices(:,2),:)];
cutNs = [n(cutIndices(:,1),:),n(cutIndices(:,2),:)];
cutKs = [K(cutIndices(:,1)),K(cutIndices(:,2))];

% Image size
imSize1 = size(S,1);

% Get the linear indices for the pixels for each cut
cutLinInds = cell(size(cuts,1),1);
for i = 1:size(cuts,1)
    inds = rayTrace(cuts(i,1:2),cuts(i,3:4));
    cutLinInds{i} = inds(:,1) + (inds(:,2)-1)*imSize1;
end

% Get the cut line segments
cutLines = reshape([cuts, nan(size(cuts,1),2)]',2,size(cuts,1)*3)';

% Find all centers that have more than one unpaired cut
unpairedCuts = cutCenters(isnan(cutCenters(:,2)),1);
[uniqueUnpairedCuts,~,ic1] = unique(unpairedCuts);
centersWithMoreThanOneUnpairedCut = uniqueUnpairedCuts(accumarray(ic1,1)>1); % Lx1

% Initialize boolian array of cuts that should be removed
cutsToRemove = false(size(cuts,1),1);

% Get a boundary offset by the inward pointing normal vectors.
offset_B = B + 0.25*n;
%     This offset boundary will be used when finding the index of an
%     intersection with the boundary. We need to offset it because it could
%     be that the boundary of a hole goes through the same pixel twice on
%     opposite sices of the contour. This would give two answers for the
%     the index of intersection at that pixel. Offseting the boundary
%     ensures that there are no repeated x-y values, and so the index of
%     intersection will be unique.

% Get the x-y centers
centers_x = centers(:,1);
centers_y = centers(:,2);



% -------------------------------------------------------------------------
% End preliminaries


% Start of computation code:
% -------------------------------------------------------------------------
% Find groups of triangles. A group of triangles are triangles that share
% an edge.

triangleGroups = triGroups(triangleList'); % cell array with elements of size of 3xLi where sum(Li) = L

if DEBUG
    Info(numel(triangleGroups)) = struct('triGroup',[],'originalCuts',[],'vertexOptimizedCuts',[],'optimizedCuts',[],'triangleCenterSearchPoints',[]);
end

% Initialize array to store any new triangle cuts
triCuts = zeros(size(triangleList,1)*3,4);
triCutCounter = 1;

for groupNum = 1:numel(triangleGroups)
    groupNum
    % Extract the current group of triangles. Remember that if this group
    % contains more than one triangle (size(triGroup,2)>1) then the
    % triangles share an edge.
    triGroup = triangleGroups{groupNum}; % 3xQ
    
    % Get a unique list of the triangle vertices
    triVrts = unique(triGroup(:))'; % 1xL
    
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
    
    % Find all cuts associated with centers of this triangle group
    % ....................................................................
    % Get paired edges
    currentCuts = sum(sum(cutCenters == permute(triVrts,[1,3,2]),3),2)==2; % Nx2xL -> Nx1
    
    % Get unpaired edges
    currentUnpairedCuts =  any(cutCenters(:,1) == triVrts,2) & isnan(cutCenters(:,2));
    
    % There should be at most one unpaired edge 
    
    % If this set of currentEdges does not contain reference to one of
    % the triangle centers, then look for unpaired edges.
    cntrsNotIncluded = ~ismember(triVrts,unique(cutCenters(currentCuts,:)));
    
    if any(cntrsNotIncluded)
        % Note that unpaired edges are always listed in the first
        % column of cutCtnrs.
        unpairedCuts = any(cutCenters(:,1) == triVrts(cntrsNotIncluded),2) & isnan(cutCenters(:,2));

        % Check to see if a center of the triangle has more than one
        % unpaired edges, and if it does, only use the edge whose midpoint
        % is closest to the mean triangle center containing the center.
        moreThanOneUnpaired = find(any(triVrts == centersWithMoreThanOneUnpairedCut,2)); %LxK -> Lx1

        for j = 1:numel(moreThanOneUnpaired)

            tmpCntr = centersWithMoreThanOneUnpairedCut(moreThanOneUnpaired(j));
            tmpCutInds = find(cutCenters(:,1) == tmpCntr);
            tmpCuts = cuts(tmpCutInds,:);
            tmpEdgesMidPoints = (tmpCuts(:,1:2) + tmpCuts(:,3:4))/2; % Lx2

            tmpTriCent = mean(triCenters(any(triGroups' == tmpCntr,2),:),1); % 1x2
            
            [~,edgeToUse] = min( sum((tmpEdgesMidPoints - tmpTriCent).^2,2) );
            tmpCutInds(edgeToUse) = [];
            unpairedCuts(tmpCutInds) = false;
        end
    else
        unpairedCuts = false(size(currentCuts));
    end

    % Combine the cuts to get the complete set.
    currentCuts = currentCuts | unpairedCuts;
    % ....................................................................
    % Section end
    
    
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
    
    newPotentialCuts = zeros(numel(sharedEdges) - sum(sharedEdges)/2,4);
    newPotentialCutInds = zeros(sum(~sharedEdges),1);
    newPotentialCutLineSegs = zeros(3*sum(~sharedEdges),2);
    
    if DEBUG
        originalCuts = zeros(sum(~sharedEdges),4);
        vertexOptimizedCuts = zeros(sum(~sharedEdges),4);
    end
    
    counter = 1;
    counter1 = 1;
    triGroup
    
    figure(1)
        clf(1)
    
    for triNum = 1:size(triGroup,2)
        triNum
        cntInds = (1:3) + (triNum-1)*3;
        isShared = sharedEdges(cntInds)>0;
        
        % Part 1:
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Create new potential cuts from the center of the triangle through
        % the mid point of any un-shared triEdges to the nearest boundary.
        
%         uv = edgeMidPoints(cntInds(~isShared),:) - triCenters(triNum,:);
%         uv = uv ./ sqrt(sum(uv.^2,2));

%         pCutEnds = triCenters(triNum,:) + MAX_CUT_LENGTH * uv(cntInds(~isShared),:);
%         
%         % Create a NaN delimited array of the new proposed cuts.
%         xTmp = [triCenters(triNum,1)*[1;1;1], pCutEnds(:,1) nan(size(pCutEnds,1),1)]';
%         yTmp = [triCenters(triNum,2)*[1;1;1], pCutEnds(:,2) nan(size(pCutEnds,1),1)]';
%         
%         % Find closest intersection of new proposed cuts with boundary.
%         [~,~,pCutIntrsctns,BIntrsctns] = intersections(xTmp(:),yTmp(:),offset_B(:,1),offset_B(:,2),0,1);
%         
%         if numel(pCutIntrsctns) > 3
%             % If there are more than three intersectinos then we need
%             % to use the intersections from the proposed cuts.
% 
%             % Find the intersection of each proposed cut that is
%             % closest to the triangle center.
%             proposedCutIntersections = ceil(pCutIntrsctns/3);
%             closestIntersection = accumarray(proposedCutIntersections, pCutIntrsctns, [3,1], @min, nan);
%             if any(isnan(closestIntersection))
%                 error('A NaN was returned when finding the closest intersection between the proposed triangle cuts and the boundary.')
%             end
%             % Compute the fractional distance along each proposed cut
%             % to the closest intersection.
%             closestIntersection = closestIntersection - floor(closestIntersection);
% 
%             % Find the xy location of the closest intersections.
%             closestIntersectionLocations = triCenters(triNum,:) + (MAX_CUT_LENGTH * closestIntersection) .* uv;
% 
%             % Find the closest boundary vertex to these intersection
%             % locations.
%             [~,pCutInds] = min(sum((offset_B - permute(closestIntersectionLocations,[3,2,1])).^2,2));
%             pCutInds = permute(pCutInds,[3,2,1]);
%         else
%             % If there are only three intersections, then we can get
%             % the intersection location directly from the boundary
%             % intersections.
%             pCutInds = round(BIntrsctns);
%         end

        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % End Part 1
        
        
        % Part 1: method 2
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
        % End Part 1: method 2
        
        numPCuts = numel(pCutInds);
        
        if DEBUG
            originalCuts(counter:counter+numPCuts-1,:) = [B(pCutInds,:), triCenters(triNum,:) .* ones(numPCuts,1)];
        end
        
        originalpCutVrts = B(pCutInds,:);
        % Part 2:
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Optimize the new potential cuts through un-shared triEdges.
        
        walls = cutIndices(~currentCuts,:);
        walls = sort([walls(:);pCutInds]);

        [pCutVrts,pCutVrtNs,~,pCutInds] = optimizeVertex(pCutInds,triCenters(triNum,:),walls,B,n,K,options);
        
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
        
        
        line(B(:,1),B(:,2),'color','w')
        line(triCenters(triNum,1),triCenters(triNum,2),'marker','x','markerSize',15,'color','r')
        line(centers_x,centers_y,'linestyle','none','marker','*','color','g')
        line(originalpCutVrts(:,1),originalpCutVrts(:,2),'linestyle','none','marker','*','color','m')
        line(pCutVrts(:,1),pCutVrts(:,2),'linestyle','none','marker','*','color','c')
%         line(triangleSearchBounds(:,1),triangleSearchBounds(:,2),'color','r','marker','o')
        drawnow;
        goDark(gcf);
        
        % Remove any search point that is outside of the search
        % triangle.
        in = inpolygon(search_r(:,1),search_r(:,2),triangleSearchBounds(:,1),triangleSearchBounds(:,2));
        search_r(~in,:) = [];
        
        line(search_r(:,1),search_r(:,2),'marker','s','linestyle','none','markersize',5,'color','y')
        
        if ~isempty(search_r)
            
            r1 = search_r - permute(pCutVrts,[3,2,1]); % Lx2xK
            d1 = sqrt(sum(r1.^2,2)); % Lx1xK
            r1 = r1 ./ d1; % Lx2xK

            if any(isShared)
                r2 = search_r - tmpCutVrts; % Lx2xR
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
        
        line(triCenters(triNum,1),triCenters(triNum,2),'marker','o','markerSize',15,'color','g')
        
        for j = 1:numPCuts
            line([pCutVrts(j,1), triCenters(triNum,1)], [pCutVrts(j,2), triCenters(triNum,2)],'Color','c','LineWidth',2)
        end
        
        if DEBUG
            Info(groupNum).triangleCenterSearchPoints = [Info(groupNum).triangleCenterSearchPoints; search_r];
        end
        
        % Save the results for later.
        triCenters(triNum,:) = round(triCenters(triNum,:));
        newPotentialCuts(counter:counter+numPCuts-1,:) = [pCutVrts, triCenters(triNum,:) .* ones(numPCuts,1)];
        newPotentialCutInds(counter:counter+numPCuts-1,:) = pCutInds;
        
        xTmp = [triCenters(triNum,1).*ones(numPCuts,1), pCutVrts(:,1) nan(numPCuts,1)]';
        yTmp = [triCenters(triNum,2).*ones(numPCuts,1), pCutVrts(:,2) nan(numPCuts,1)]';
        
        newPotentialCutLineSegs(counter1:counter1+numel(xTmp)-1,:) = [xTmp(:), yTmp(:)];
        
        counter = counter + numPCuts;
        counter1 = counter1 + numel(xTmp);
    end % for each triangle
    
    if DEBUG
        Info(groupNum).originalCuts = originalCuts;
        Info(groupNum).vertexOptimizedCuts = vertexOptimizedCuts;
    end
    
    % ....................................................................
    % * Determine if any of the new potential cuts intersect with any cut
    %   not associated with this triangle group, and if they do, do not
    %   consider this triangle group.
    
    % Get the intersections between the new proposed edges and the old
    % set of edges.
    [~,~,intrsctns] = intersections(cutLines(:,1),cutLines(:,2),newPotentialCutLineSegs(:,1),newPotentialCutLineSegs(:,2));
    intersectedEdgeInds = unique(ceil(intrsctns/3));
    intersectedEdgeInds
    currentCuts
    cuts
    cutCenters
    find(currentCuts)
    line(cutLines(:,1),cutLines(:,2),'color','m')
    line(newPotentialCutLineSegs(:,1),newPotentialCutLineSegs(:,2),'color','r')
    disp here
    if any(~ismember(intersectedEdgeInds, find(currentCuts)))
        continue;
    end
    % ....................................................................

    % ....................................................................
    % Get the average new potential cut distance, the averge dot product of
    % the vertex' normals with the cut direction, and the average curvature
    % of the cut vertices.
    r1 = newPotentialCuts(1:counter-1,3:4) - newPotentialCuts(1:counter-1,1:2); % Lx2
    n1 = n(newPotentialCutInds(1:counter-1),:);
    k1 = K(newPotentialCutInds(1:counter-1));
    
    d1 = sqrt(sum(r1.^2,2)); % Lx1
    r1 = r1 ./ d1;
    a1 = sum(r1 .* n1,2);
    
    newD = mean(d1);
    newA = mean(a1);
    newK = mean(k1);
    % ....................................................................
    
    % ....................................................................
    % Get the average cut distance, the average dot product of the vertex'
    % normals with the cut directions, and the average curvatrue of the cut
    % vertices.
    r1 = cuts(currentCuts,1:2);
    n1 = cutNs(currentCuts,1:2);
    
    r2 = cuts(currentCuts,3:4);
    n2 = cutNs(currentCuts,3:4);
    
    k = cutKs(currentCuts,:);
    
    r = r1-r2;
    d = sqrt(sum(r.^2,2));
    r = r ./ d;
    
    a1 = sum(-r .* n1,2);
    a2 = sum( r .* n2,2);
    
    oldD = mean(d);
    oldA = mean([a1;a2]);
    oldK = mean(k(:));
    
    % ....................................................................
    
    % ....................................................................
    % * Create new potenial cuts between each of the triangle centers that
    %   share an edge

    for sharedEdgeNum = 1:max(sharedEdges)
       pairedCenters = ceil(find(sharedEdges==sharedEdgeNum)/3);
       newPotentialCuts(counter,:) = [triCenters(pairedCenters(1),:), triCenters(pairedCenters(2),:)];
       counter = counter + 1;
       % Create a cut between the two pairedCenters.
    end
    % ....................................................................
    
    if DEBUG
        Info(groupNum).optimizedCuts = newPotentialCuts;
    end
    
    % ....................................................................
    % * Overlap the new potential cuts with the edge image. Overlap the
    %   cuts associated with this triangle group with the edge image. Keep
    %   the set of edges (new potential cuts or cuts associated with this
    %   triangle group) with the larges mean overlap
    
    newPotentialCutLinInds = cell(size(newPotentialCuts,1),1);
    for cutNum = 1:size(newPotentialCuts,1)
        inds = rayTrace(newPotentialCuts(cutNum,1:2),newPotentialCuts(cutNum,3:4));
        newPotentialCutLinInds{cutNum} = inds(:,1) + (inds(:,2)-1)*imSize1;
    end

    newO = mean(S(cat(1,newPotentialCutLinInds{:})));
    oldO = mean(S(cat(1,cutLinInds{currentCuts})));
    
    % Normalized old and new scores (by type) by their mean
    scores = [newO,oldO;newA,oldA;newK,oldK;newD,oldD]
    scores = scores ./ mean(scores,2)
    score = sum(scores(1:3,:),1)% ./ scores(4,:);
%     newScore = (newA + newK + newO)/newD;
%     oldScore = (oldA + oldK + oldO)/oldD;
%     [newScore,oldScore]
    if score(1) > score(2)
        cutsToRemove(currentCuts) = true;
        triCuts(triCutCounter:triCutCounter+size(newPotentialCuts,1)-1,:) = newPotentialCuts;
        triCutCounter = triCutCounter + size(newPotentialCuts,1);
    end
    
end


end



function triangleList = nonBoundaryIntersectingTriangles(B,centers)
% NONBOUNDARYINTERSECTINGTRIANGLES  Compute a list of traingles with
% vertices given by CENTERS whose edges do not overlap with the boundary B.

% Get the daulaunay triangulation
triangleList = delaunay(centers);

% Get the indices for each edge of the triangulation
triEdges = [triangleList(:,[1,2]); triangleList(:,[2,3]); triangleList(:,[3,1])];
triEdges = unique(sort(triEdges, 2, 'ascend'),'rows');

% From the line segments of each edge as a NaN delimited list.
xVerts = [centers(triEdges(:,1),1) centers(triEdges(:,2),1) nan(size(triEdges,1),1)]';
yVerts = [centers(triEdges(:,1),2) centers(triEdges(:,2),2) nan(size(triEdges,1),1)]';

% Get the intersections of the triangulation edges with the boundary.
[~,~,intrsctns] = intersections(xVerts(:),yVerts(:),B(:,1),B(:,2),0,1);

% Get the indices of each edge that intersects with the boundary.
badEdges = triEdges(unique(ceil(intrsctns/3)),:);

% Remove all triangles with an edge that intersects the boundary.
for i = 1:size(badEdges,1)
    triangleList(any(triangleList == badEdges(i,1),2) & any(triangleList == badEdges(i,2),2),:) = [];
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


% Getting the search range around each unpairedCut vertex requires some
% care. We cannot simply use B(i-searchRange:i+searchRange), but the search
% range 1) wrap around each contour, and 2) only contain indices of the
% same contour as the upairedCut vertex. Furthermore, the searchRange
% should not cross any other cut vertex.

% Get search range.
SEARCH_RADIUS = options.searchRadius;
WIGNER_SEITZ_RADIUS = options.WIGNER_SEITZ_RADIUS;
maxK = 1/WIGNER_SEITZ_RADIUS;


% Put limits on the curvature;
K(K>maxK) = maxK;
K(K<-maxK) = -maxK;


% Get the contours of the boundary
nanLocs = isnan(B(:,1));
BcontourNumber = cumsum(nanLocs) + 1;
Binds = (1:size(B,1))';
% 
% cutIndices = cutIndices(unpairedCuts,:);
% cutIdxContourNumbers = reshape(BcontourNumber(cutIndices),size(cutIndices));


vertexContourNumbers = BcontourNumber(vertexIndices);
wallContourNumbers = BcontourNumber(wallIndices);


walls = accumarray(wallContourNumbers(:), wallIndices(:),[sum(nanLocs)+1,1],@(x) {sort(x)},[]);
contourInds = accumarray(BcontourNumber(~nanLocs),Binds(~nanLocs),[sum(nanLocs)+1,1],@(x) {x},[]);
% walls
optimizedVertexIndices = zeros(size(vertexIndices,1),1);

for i = 1:size(vertexIndices,1)

    % Get the search indices for the two vertices of the current cut
    c1 = vertexContourNumbers(i);
    [searchInds1, ind1] = getSearchInds(vertexIndices(i),contourInds{c1},walls{c1},SEARCH_RADIUS);

%     fprintf('-----------------------------------------------------------\n');
%     fprintf('%d\n',i)
%     walls{1}
%     walls{2}
%     unpairedCutsIdx(i,:)
%     searchInds1(:)'
%     searchInds2(:)'
    % Get the locations and normals over each sarch range.
    r1 = B(searchInds1,:); % N1 x 2
    n1 = n(searchInds1,:);
    k1 = K(searchInds1,:);
% [searchInds1(:),k1]

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
    
%     [n1_ind,n2_ind,searchInds1(n1_ind),searchInds2(n2_ind),ind1,ind2]
    
    % Update the wall positions.
    walls{c1}(ind1) = searchInds1(ind);

    % Save the new indices
    optimizedVertexIndices(i) = searchInds1(ind);
    
end

cut = B(optimizedVertexIndices,:);
cutN = n(optimizedVertexIndices,:);
cutK = K(optimizedVertexIndices);

end