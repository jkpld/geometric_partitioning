function [cut,cutN,cutK,optimizedCutIndices] = optimizeCuts(cutIndices,B,n,K,options)
% OPTIMIZECUTS  Search in the local region of the vertices that make up a
% cut for new vertices that minimize the cut distance, maximize the dot
% product between the vertex normals and the cut direction, and maximize
% the curvature.
%
% Input parameters:
%
% cutIndices - Nx2 array where cutIndices(i,1) and cutIndices(i,2) give the
%              indices into B for the the i'th cut' vertices.
%
% B - Mx2 array with the boundary
%
% n - Mx2 array with the boundary unit normals (pointing inward)
%
% K - Mx1 array with boundary curvature
%
% options - structure with one required field
%           Search_Radius - The radius of the search range over which the
%           unpaired cut will be optimized.
%
% Output parameters:
%
% cut - Nx4 array where cut(i,1:2) and cut(i,3:4) give the vertices of the 
%       i'th cut. 
%
% cutN - Nx4 array where cutN(i,1:2) and cutN(i,3:4) give the unit normals
%        (pointing inward) of the i'th cut vertices.
%
% cutK - Nx2 array where cutK(i,1) and cutK(i,2) give the boundary
%        curvature at the i'th cut vertices.
%
% optimizedCutIndices - Nx2 array where optimizedCutIndices(i,:) gives the
%                       indices of the vertices of the i'th cut in the
%                       boundary.
%
% A cut is optimized by searching for a new vertices within searchRadius of
% the current vertices such that the distance of the cut is minimized, the
% dot product of the cut' vertex normals with the cut direction is
% maximized, and the curvature is maximized. The search region is not
% allowed to cross any other cut vertex.
%
% See also GETSEARCHINDS CREATEOBJECTPARTITIONS

% James Kapaldo
% 2016-10-06


% Getting the search range around each unpairedCut vertex requires some
% care. We cannot simply use B(i-searchRange:i+searchRange), but the search
% range 1) wrap around each contour, and 2) only contain indices of the
% same contour as the upairedCut vertex. Furthermore, the searchRange
% should not cross any other cut vertex.

% Get search range.
SEARCH_RADIUS = options.Search_Radius;


% Put limits on the curvature - Not used anymore
% WIGNER_SEITZ_RADIUS = options.Wigner_Seitz_Radius;
% maxK = 1/WIGNER_SEITZ_RADIUS;
% K(K>maxK) = maxK; 
% K(K<-maxK) = -maxK;
K(K<0) = K(K<0)*5;

% Get the contours of the boundary
nanLocs = isnan(B(:,1));
BcontourNumber = cumsum(nanLocs) + 1;
Binds = (1:size(B,1))';

cutIdxContourNumbers = BcontourNumber(cutIndices);

if size(cutIndices,1)==1
    cutIdxContourNumbers = cutIdxContourNumbers';
end

walls = accumarray(cutIdxContourNumbers(:), cutIndices(:),[sum(nanLocs)+1,1],@(x) {sort(x)},[]);
contourInds = accumarray(BcontourNumber(~nanLocs),Binds(~nanLocs),[sum(nanLocs)+1,1],@(x) {x},[]);

optimizedCutIndices = zeros(size(cutIndices,1),2);

for i = 1:size(cutIndices,1)
    
    % Get the search indices for the two vertices of the current cut
    c1 = cutIdxContourNumbers(i,1);
    c2 = cutIdxContourNumbers(i,2);
    [searchInds1, ind1] = getSearchInds(cutIndices(i,1),contourInds{c1},walls{c1},SEARCH_RADIUS);
    [searchInds2, ind2] = getSearchInds(cutIndices(i,2),contourInds{c2},walls{c2},SEARCH_RADIUS);

    if numel(ind1) ~= 1 || numel(ind2) ~= 1
%         cutIndices(:)
%         for j = 1:numel(walls)
%             walls{j}
%         end
        error('optimizeCuts:duplicateWalls', 'There are at least two indices that are the same in the WALLS array.')
    end
    
    % Get the locations and normals over each sarch range.
    r1 = B(searchInds1,:); % N1 x 2
    n1 = n(searchInds1,:);
    k1 = K(searchInds1,:);

    r2 = B(searchInds2,:); % N2 x 2
    n2 = n(searchInds2,:);
    k2 = K(searchInds2,:);

    % Get the vector between each possible search location
    r = r1 - permute(r2,[3,2,1]); % N1 x 2 x N2
    d = sqrt(sum(r.^2,2)); % N1 x 1 x N2
    r = r ./ d; % N1 x 2 x N2
    
    % Get the dot products between the vertex normals and the cut
    % directions
    a1 = sum( -r .* n1, 2); % N1 x 1 x N2
    a2 = sum( r .* permute(n2,[3,2,1]),2); % N1 x 1 x N2
    
    % Get curvature combinations
    k = k1 + permute(k2,[3,2,1]); % N1 x 1 x N2
    
    % Optimize
    valueToOptimize = permute((a1+a2+k)./d,[1,3,2]); % N1 x N2
    [~,ind] = max(valueToOptimize(:));
    
    n1_ind = rem(ind-1, numel(searchInds1)) + 1; % N1 ind
    n2_ind = (ind - n1_ind)/numel(searchInds1) + 1; % N2 ind
    
    % Update the wall positions.
    walls{c1}(ind1) = searchInds1(n1_ind);
    walls{c2}(ind2) = searchInds2(n2_ind);
    
    walls{c1} = sort(walls{c1});
    if c1 ~= c2
        walls{c2} = sort(walls{c2});
    end
        
    

    % Save the new indices
    optimizedCutIndices(i,:) = [searchInds1(n1_ind),searchInds2(n2_ind)];
    
end

optimizedCutIndices = sort(optimizedCutIndices,2);
cut = [B(optimizedCutIndices(:,1),:),B(optimizedCutIndices(:,2),:)];
cutN = [n(optimizedCutIndices(:,1),:),n(optimizedCutIndices(:,2),:)];
cutK = [K(optimizedCutIndices(:,1)),K(optimizedCutIndices(:,2))];

end
% optimizeCuts : changeLog
% 2016-10-15 : Re-written. Old version commented out below.
% 2016-10-16 : Added in curvature to the optimization. Renamed to
%              optimizeCuts. Added ouput for the cuts, cut' normals, and
%              cut' curvatures. Removed unpairedIndices from input.
% 2016-10-20 : removed limit on max curvature
% 2016-10-21 : removed limit on min curvature and infact increased the
%              neagive effect of negative curvatures my multiplying by 5.
%              this should help with the problem that a cut slides to a
%              corner of a nucleous for the benafit of being shorter.

% function [optimizedCuts,optimizedCutVertNorms,optimizedCutVertCurvatures,cutVertices] = optimizeUnpairedCuts(unpairedCuts,B,n,curvature,cutVertices,options)
% % optimizeUnpairedCuts will shift the vertices of each unpaired cut to
% % maximize the dot product between the boundary normals at the vertex and
% % the cut direction and minimize the distance between the two cut vertices.
% %
% % Input parameters:
% %
% % unpairedCuts - Nx4 array where N is the number of unpaired cuts to
% %                optimize. unpairedCuts(i,1:2) will give the *indices* of
% %                the two vertices that make up the i'th cut.
% %                unpairedCuts(i,3:4) give the contour number of the
% %                vertices that make up the i'th cut.
% %
% % B - A cell array of the boundary contours, as returned by
% %     splitContourAndGetIndices.
% %
% % n - A cell array of boundary contour normals, as returened by
% %     splitContourAndGetIndices.
% %
% % curvature - A cell array of boundary contour curvatures, as returned by
% %             splitContourAndGetIndices.
% %
% % cutVertices - A cell array with the indices of all of the cut vertices,
% %               as returned by splitContourAndGetIndices.
% %
% % options - structure with possible fields
% %             searchRadius - The radius of the search range for optimizing
% %                            the cuts in units of boundary pixels.
% %
% % Output parameters:
% %
% % optimizedCuts - Nx4 array where N is the number of unpaired cuts.
% %                 optimizedEdges(i,1:2) and optimizedEdges(i,3:4) give the
% %                 *x-y* locations of the two vertices that make up the
% %                 optimized cut.
% %
% % optimizedCutsVertNorms - Nx4 array where N is the number of unpaired 
% %                          cuts. optimizedCutsVertNorms(i,1:2) and
% %                          optimizedCutsVertNorms(i,3:4) give the normal 
% %                          vectors at the each of the two vertices that 
% %                          make up the optimized cut.
% %
% % optimizedCutVertCurvatures - Nx2 array where N is the number of unapired
% %                              cuts. optimizedCutsVertCurvatures(i,1) and
% %                              optimizedCutsVertCurvatures(i,2) give the
% %                              curvatures of the two vertices that make up
% %                              the optimized cut.
% %
% % cutVertices - The same as the input cutVertices; however, the indices of
% %               the vertices that changed have been updated to the
% %               optimized indices.
% %
% % See also SPLITCONTOURANDGETINDICES, GETSEARCHINDS, CREATEOBJECTPARTITIONS
% 
% % James Kapaldo
% % 2016-10-06
% 
% 
% % Optimize unpaired cuts by moving their end points within 10 pixels of
% % their current possition and both maximized the dot product between the
% % boundry normals and the cut and minimizing the length of the cut.
% %
% % Note that the range of pixels each cut end point can move to will be
% % further bouned by existing paired edges. (We dont want to move an
% % unpaired edge accross another cut.)
% 
% % Get search range.
% SEARCH_RADIUS = options.searchRadius;
% 
% % Get the length of each contour.
% M = cellfun(@length,B);
% 
% % Get the number of unpairedCuts
% N = size(unpairedCuts,1);
% 
% % Initialize arrays
% optimizedCuts = zeros(N,4);
% optimizedCutVertNorms = zeros(N,4);
% optimizedCutVertCurvatures = zeros(N,2);
% 
% for i = 1:N
% 
%     % Get the search indices for each vertex
%     [searchInds1,ind1] = getSearchInds(unpairedCuts(i,1),M(unpairedCuts(i,3)),cutVertices{unpairedCuts(i,3)},SEARCH_RADIUS);
%     [searchInds2,ind2] = getSearchInds(unpairedCuts(i,2),M(unpairedCuts(i,4)),cutVertices{unpairedCuts(i,4)},SEARCH_RADIUS);
% 
%     r1 = B{unpairedCuts(i,3)}(searchInds1,:);
%     n1 = n{unpairedCuts(i,3)}(searchInds1,:);
%     k1 = curvature{unpairedCuts(i,3)}(searchInds1);
%     
%     r2 = B{unpairedCuts(i,4)}(searchInds2,:);
%     n2 = n{unpairedCuts(i,4)}(searchInds2,:);
%     k2 = curvature{unpairedCuts(i,4)}(searchInds2);
% 
%     r = r1 - permute(r2,[3,2,1]);
%     d = sqrt(sum(r.^2,2));
%     r = r ./ d; 
%     
%     a1 = sum( -r .* n1, 2);
%     a2 = sum( r .* permute(n2,[3,2,1]),2);
%     
%     valueToOptimize = permute((a1+a2)./d,[1,3,2]);
%     [~,ind] = max(valueToOptimize(:));
%     
%     j = rem(ind-1, numel(searchInds1)) + 1;
%     k = (ind - j)/numel(searchInds1) + 1;
% 
%     optimizedCuts(i,:) = [r1(j,:), r2(k,:)];
%     optimizedCutVertNorms(i,:) = [n1(j,:), n2(k,:)];
%     optimizedCutVertCurvatures(i,:) = [k1(j), k2(k)];
%     
%     % Update the cut vertex indices.
%     cutVertices{unpairedCuts(i,3)}(ind1) = searchInds1(j);
%     cutVertices{unpairedCuts(i,4)}(ind2) = searchInds2(k);
%     
%     % At this point it might also be good to fix
%     % partition{unPairedPartitionNumber(i)} to that the new cut locations
%     % is taken into account. However, this will be skipped for now as the
%     % partitions may not be needed. 
% end
% 
% end
