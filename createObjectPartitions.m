function [cuts,ncuts,cutCntrs,cutIndices,varargout] = createObjectPartitions(B,n,centers,options)
% CREATEOBJECTPARTITIONS  Partition a boundary into disjoint regions.
% The number of partitions will be given by the number of CENTERS.
%
% Input parameters:
%
% B - M x 2 array giving the boundary coordinates.
%
% n - M x 2 array giving the boundary unit normals (pointing inward).
%
% centers - N x 2 array giving the centers each partition should be built
%           arround.
%
% options - structure containing additional options.
%           Required fields:
%               Max_Radius - This gives the maximum allowable distance
%                            between a boundary and a center.
%               Min_Angle - The minimum *dot product* between the boundary
%                           element normal (point inwards) and the line
%                           connecting the boundary element to a center.
%               Use_GPU - logical, if true, then a GPU will be used to
%                         speed up intersection calculations.
%               Debug - logical, if true, then the optional Info structure
%                       will be output. If false, it will not be output.
%           Optional fields: none
%
% Output parameters:
%
% cuts - K x 4 array giving the K cuts that break up the boundary into the
%        partitions. K(i,1:2) gives one vertex of a cut and K(i,3:4) give
%        the other vertex of a cut.
%
% ncuts - K x 4 array giving the normal vectors at each cut vertex.
%
% cutCntrs - K x 2 array giving the center(s) that each cut is assigned to.
%            If a cut k splits a region between center i and center j, then
%            cutCntrs(k,:) = [i,j]. If one side of the edge does not have a
%            center, then cutCntrs(k,:) = [i, NaN].
%
%
% Optional output
% Info - a structure contaning information potentially usefull for
%        debuging.
%           assignedCenter - Nx1 array with the assigned center of each
%           boundary element.
%
%           Number_Of_Edge_Optimization_Iterations - the number of times
%           the while loop that removes crossing boundary-center
%           assignments is run.
%
%           partitions - 1 x N cell array containing B_i x 2 arrays with
%           the boundaries for each partition.
%
% Notes: 
%
% * This function assumes the boundary came from an image. This means that
% all of the x,y values in B are all integers with adjacent boundary
% verticies being at most sqrt(2) apart from each other.
% * Holes may be included in the boundary. They should be delimeted by a row
% of nan's. (The array of boundary tangents must also have the tangents for
% the hole with a row of nan delimeters.)
% * The outer boundary should be clockwise oriented. Holes should be
% counter-clockwise oriented. --- These conditions will not be checked, and
% if they are not met, then unexpected results could occure.
% * If Max_Radius is very large, then the results should still correctly
% calculate the partitions; however, the time of computation will increase.
% Thus, make Max_Radius large enough for the largest single object, but not
% larger.
%
% See also CALCULATECUTS

% James Kapaldo
% 2016-10-05

% Output Info if requrested.
debug = options.Debug;

MAX_RADIUS = options.Max_Radius;
MIN_ANGLE = options.Min_Angle;
useGPU = options.Use_GPU;
SMALL_REGION_SIZE = 5;

M = size(B,1);
N = size(centers,1);

% Get the vectors between each boundary vertex and each center. Compute
% the lengt and direction of these vectors.
dr = centers - permute(B,[3,2,1]); % Nx2xM
r_length = sqrt(sum(dr.^2,2)); % Nx1xM
dr_uv = dr ./ r_length; % Nx2xM

% Get the end points of the lines that go from each boundary element to
% each center location. These lines will be offset from the boundary or the
% center by a small amount so that they did not intersect with the
% boundary, and good lines to not intersect with each other.
tmpcnts = centers - 0.2*dr_uv; % Nx2xM -- start points
tmpbndry = permute(B,[3,2,1]) + 0.2*dr_uv; % Nx2xM -- end points

% Align pages into N*Mx2
tmpcnts = reshape(permute(tmpcnts,[2,1,3]),2,N*M)';
tmpbndry = reshape(permute(tmpbndry,[2,1,3]),2,N*M)';

% Get the dot product (angle) between the the boundary normals and the
% vectors to the centers.

drdbndry = sum(dr_uv .* permute(n,[3,2,1]),2); % Nx1xM

% Require the distance between the boundary and the centers to be less than
% maxRadius and require that the dot product be greater than 0.

validDist = r_length < MAX_RADIUS; % Nx1xM
validAngle = drdbndry > MIN_ANGLE; % Nx1xM
validPoints = validDist & validAngle; % Nx1xM

validPointsInd = find(validPoints(:)); % (<=N*M)x1

if (numel(validPointsInd) < 0.5*N)
    error('createObjectPartitions:tooFewValidPoints','The number of valid points is less than half the number of boundary vertices. Try increasing Max_Radius and/or decreasing Min_Angle.')
end

lineSegments = nan(3*numel(validPointsInd),2);
lineSegments(1:3:end,:) = tmpcnts(validPointsInd,:);
lineSegments(2:3:end,:) = tmpbndry(validPointsInd,:);

% Find intersections between line segments and boundary edges.
[~,~,intsctns] = intersections(lineSegments(:,1),lineSegments(:,2),B(:,1),B(:,2),0,0,useGPU);

% Set all points that intersect with the boundary to be not valid.
validPoints(validPointsInd(unique(ceil(intsctns/3)))) = false;
invalidPoints = permute(all(~validPoints,1),[3,2,1]); % Mx1

% Set invalid points to have an angle of 0 and a length of inf.
drdbndry(~validPoints) = 0;
r_length(~validPoints) = inf;

% Assign a boundry vertex to the center that has the largest
% drdbndry/r_length. This should both maximize the the dot product and
% minimize the distance between the boundary and center.

[selectedScores,selectedCenter] = max(drdbndry./r_length,[],1);
selectedCenter = permute(selectedCenter,[3,2,1]); % Mx1
selectedScores = permute(selectedScores,[3,2,1]); % Mx1

% mophologically close the set of invalid points. This will remove small
% segments of boundary pixels that have been assigned a center but are
% surrounded by invalid points.
% invalidPoints = imclose(invalidPoints,ones(10,1));

% Replace invalid points with nan
selectedCenter(invalidPoints) = nan; % Mx1
selectedScores(invalidPoints) = nan; % Mx1

if debug
    Info.initialAssignedCenter = selectedCenter;
end

counter = 1;

while 1
    % Now we need to intersect each winning line from the boundary to the
    % center with each other. These line are not allowed to intersect.
    
%     fprintf('counter = %d\n',counter)

    % Get the linear indicies for the selectedCenter x,y positions in
    % tmpcnts and tmpbndry.
    inds = (selectedCenter + (0:M-1)'*N) + [0, M*N];
    
    % Remove invalid center assignements.
    valid = isfinite(inds(:,1));
    inds(~valid,:) = [];
    scores = selectedScores(valid);
    valid = find(valid);
        
    % Turn the assignments into nan delimeted line segments for
    % intersection calculation.
    lineSegments = nan(3*size(inds,1),2);
    lineSegments(1:3:end,:) = tmpcnts(inds);
    lineSegments(2:3:end,:) = tmpbndry(inds);
    
    % Find the intersections.
    [~,~,i1,i2] = intersections(lineSegments(:,1),lineSegments(:,2),lineSegments(:,1),lineSegments(:,2),0,0,useGPU);

    % If there are no intersections, then we are done.
    if isempty(i1)
        break;
    end
    
    % Find the paired intersecting edge indices.
    iE = [ceil(i1/3),ceil(i2/3)]; % index of the edges that intersect
    
    % Give each edge in an intersection a score. Consider edge e1
    % intersecting with edge e2. In iE, there will be two rows, [e1, e2]
    % and [e2, e1]. We want to remove one of these edges that results in
    % the "best" set of edges. The scores previously calculated were
    % maximized for good bndry-center assignements. Use these scores to
    % calculate a new score for each row of iE. The score will be the ratio
    % of the first column score to the second column score.
    %
    % S([e1,e2]) = S(e1)/S(e2);
    % S([e2,e1]) = S(e2)/S(e1);
    %
    % We will then accumulate the first column of iE to give a single
    % unique array with each edge listed only once, and we will have a
    % single score array with the accumulated scores for each edge. We will
    % then iteratively remove the edge with the lowest score until there
    % are no more intersections.
    
    Score = scores(iE(:,1))./scores(iE(:,2));
    
    iEinds = unique(iE(:,1));
    binEdges = [iEinds(1)-0.5; ...
                iEinds(1:end-1) + diff(iEinds)/2; ...
                iEinds(end)+0.5];
    [nI,bin] = histcountsmex(iE(:,1),binEdges);
    totalScore = accumarray(bin,Score,size(iEinds),@mean,0);

    % We know need to remove edges until there are no intersections left.
    % There are two methods we could use. Iteratively remove the edge with
    % the most intersections, or iteratively remove the longest edge with
    % intersections.
    edgesToRemove = zeros(numel(iEinds)-1,1);
    counter2 = 1;

    % ********************************************************************
    % Method 3 ------------------------------------------------- the best
    %
    % Iteratively remove the edge with the lowest score
    while any(nI >= 1)
        [~,mostIntsInd] = min(totalScore);
        mostIntersections = iEinds(mostIntsInd);
        
        edgesToRemove(counter2) = mostIntersections;
        
        iE(any(iE==mostIntersections,2),:) = [];
        Score = scores(iE(:,1))./scores(iE(:,2));
        [nI,bin] = histcountsmex(iE(:,1),binEdges);

        totalScore = accumarray(bin,Score,size(iEinds),@mean,inf);
        
        counter2 = counter2 + 1;
    end
    edgesToRemove(counter2:end) = [];
    % End Method 3 ------------------------------------------------------- 
    
    % Quit if there are no edges to remove.
    if isempty(edgesToRemove)
        break;
    end
    
    % Now set these edges to invalid and recompute
    validPoints((valid(edgesToRemove)-1)*N + selectedCenter(valid(edgesToRemove))) = false;
    
    invalidPoints = permute(all(~validPoints,1),[3,2,1]); % Mx1
    
    % Set invalid points to have an angle of 0 and a length of inf.
    drdbndry(~validPoints) = 0;
    r_length(~validPoints) = inf;
    
    % Assign a boundry vertex to the center that has the largest
    % drdbndry/r_length. This should both maximize the the dot product and
    % minimize the distance between the boundary and center.
    
    [~,selectedCenter] = max(drdbndry./r_length,[],1);
    selectedCenter = permute(selectedCenter,[3,2,1]); % Mx1
    
    % mophologically close the set of invalid points. This will remove small
    % segments of boundary pixels that have been assigned a center but are
    % surrounded by invalid points.
%     invalidPoints = imclose(invalidPoints,ones(10,1));
    
    % Replace invalid points with nan
    selectedCenter(invalidPoints) = nan; % Mx1
    counter = counter + 1;
end

% -----------------------------------------------------------------------
% Replace small regions of invalid points surrounded by the same
% selectedCenter with the selectedCenter. 
% This will fix problems when the assignment line from a boundary element
% to a center intersects with the boundary elements just adjacent to it.

nanLocs = [0;find(isnan(B(:,1)));size(B,1)+1];

for i = 1:numel(nanLocs)-1
    % Get the current contour
    regionInds = nanLocs(i)+1:nanLocs(i+1)-1;
    tmpInvalidPoints = invalidPoints(regionInds);
    tmpSelectedCenter = selectedCenter(regionInds);
    tmpM = numel(regionInds);

    % If there are invalid points in this contour then iterate over
    % each connected group.
    if any(tmpInvalidPoints)

        % Get the groups
        edgs = find( tmpInvalidPoints - circshift(tmpInvalidPoints,-1,1));
        edgs_p1 = circshift(edgs,-1,1);
        edgs_m1 = circshift(edgs,1,1);

        if tmpInvalidPoints(1)
            % If the first point is invalid then there is the special
            % case where a region of invalid points stradles the
            % boundary start/end.
            if ((edgs(1)+M-edgs_m1(1)) < SMALL_REGION_SIZE) && (tmpSelectedCenter(edgs(1)+1)==tmpSelectedCenter(edgs_m1(1)))
                tmpSelectedCenter(edgs_m1(1)+1:end) = tmpSelectedCenter(edgs_m1(1));
                tmpSelectedCenter(1:edgs(1)) = tmpSelectedCenter(edgs(1)+1);
            end
            % Go over all other invalid regions
            for j = 2:2:numel(edgs)-2
                if (edgs_p1(j)-edgs(j) < SMALL_REGION_SIZE) && (tmpSelectedCenter(edgs(j))==tmpSelectedCenter(mod(edgs_p1(j),tmpM)+1))
                    tmpSelectedCenter(edgs(j)+1:edgs_p1(j)) = tmpSelectedCenter(edgs(j));
                end
            end
        else
            % No invalid region stradles the boundary start/end so just
            % go over the regions.
            for j = 1:2:numel(edgs)-1
                if (edgs_p1(j)-edgs(j) < SMALL_REGION_SIZE) && (tmpSelectedCenter(edgs(j))==tmpSelectedCenter(mod(edgs_p1(j),tmpM)+1))
                    tmpSelectedCenter(edgs(j)+1:edgs_p1(j)) = tmpSelectedCenter(edgs(j));
                end
            end
        end

        selectedCenter(regionInds) = tmpSelectedCenter;
    end
end

% Compute cuts
if debug
    [cuts,ncuts,cutCntrs,cutIndices,partitions,CalculateCutsInfo] = calculateCuts(B,n,selectedCenter,centers);
else
    [cuts,ncuts,cutCntrs,cutIndices,partitions] = calculateCuts(B,n,selectedCenter,centers);
end


if debug
    Info.assignedCenter = selectedCenter;
    Info.Number_Of_Edge_Optimization_Iterations = counter-1;

    CalculateCutsInfo.partitions = partitions;
    Info.CalculateCuts = CalculateCutsInfo;
    varargout{1} = Info;
end

end
% createObjectPartitions : changeLog
% 2016-10-14: Added in code to remove small regions of invalid points
%             surrounded by the same selected center. 
% 2016-10-21: Fixed but in code to remove the small regions of invalid
%             points.
% 2016-10-25: reduced small region size to 5.
% 2016-10-28: moved partitions to info. renamed required fields in options













% ********************************************************************
% Method 1 ----------------------------------------------- works okay
%
% Iteratively remove the edge with the most intersections until there
% are no more intersections.
% while any(nI>1)
%     [~,mostIntsInd] = max(nI);
% 
%     mostIntersections = iEinds(mostIntsInd);
% 
%     edgesToRemove(counter2) = mostIntersections;
% 
%     iE(any(iE==mostIntersections,2),:) = [];
%     nI = histcountsmex(iE(:,1),binEdges);
% 
%     counter2 = counter2 + 1;
% end
% edgesToRemove(counter2:end) = [];
% End Method 1 -------------------------------------------------------
    
% ********************************************************************
% Method 2 --------------------------------------- does not work well
%
% Iteratively remove the longest edge that has at least 2
% intersections.
% haveIntersections = numIntersections>1;
% while any(haveIntersections)
%     intrsctnInds = intersectingEdgeInds(haveIntersections);
%     rc = tmpcnts(valid(intrsctnInds),:);
%     rb = tmpbndry(valid(intrsctnInds),:);
% 
%     d = sum((rc-rb).^2,2);
%     [~,maxInd] = max(d);
%     edgesToRemove(counter2) = intrsctnInds(maxInd);
% 
%     intersectingEdges(any(intersectingEdges==intrsctnInds(maxInd),2),:) = [];
%     numIntersections = histcountsmex(intersectingEdges(:,1),binEdges);
%     haveIntersections = numIntersections>1;
% 
%     counter2 = counter2 + 1;
% end
% edgesToRemove(counter2:end) = [];
% End Method 2 -------------------------------------------------------