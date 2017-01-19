function [cutA, vrtA] = computeCutAdjacency(cutInds,B,adjacencyDistance)
% COMPUTECUTADJACENCY  Determine which cuts are adjacent to each
% other and return the adjacency matrix.
%
% Two cuts are adjacent if they have vertices within ADJACENCYDISTANCE of
% each other
%
% cutA - the cut adjacnecy matrix with size M x M, where M is the number of
% cuts, M = size(cutInds,1);
%
% vrtA - the vertex adjacency matrix of size 2*M x 2*M, the numbering is
% for the reshaped cuts reshape(cutInds',2*M,1), where M = size(cutInds,1)
%
% If ADJACENCYDISTANCE has more than one element, then cutA and vrtA will
% both be cell arrays of size 1 x numel(ADJACENCYDISTANCE). The ith element
% of the cell arrays will hold the cut/vertex adjacency matrix for the ith
% ADJACENCYDISTANCE.

% James Kapaldo
% 2016-10-20


if size(cutInds,1) == 1 
    cutA = 0;
    vrtA = zeros(2);
    return;
end

ADJACENCY_DISTANCE = adjacencyDistance;


M = size(cutInds,1);
cutInds = reshape(cutInds',2*M,1);


% Get the contours of the boundary
nanLocs = isnan(B(:,1));

if any(nanLocs)
    
    % If there are multiple contours then we need offset the indices of the
    % input array so that that each contour starts at 1.
    
    nanLocs = double(nanLocs);
    BcontourNumber = cumsum(nanLocs) + 1;

    % Lengths of each contour
    contourLengths = [find(nanLocs,1,'first');diff([find(nanLocs);size(B,1)])];

    % Get indice offsets
    nanLocs(logical(nanLocs)) = contourLengths(1:end-1);
    offsets = cumsum(nanLocs);

    % Get the contour number of each index
    triCent = isnan(cutInds);
    cutContour = nan(size(cutInds));
    cutContour(~triCent) = BcontourNumber(cutInds(~triCent));
    
    % Let the index of each contour start at 1
    cutInds(~triCent) = cutInds(~triCent) - offsets(cutInds(~triCent));
    
    % Create an array the same size as cutInds with the contour length
    contourLength = nan(size(cutContour));
    contourLength(~triCent) = contourLengths(cutContour(~triCent)) - 1;
    
else
    contourLength = size(B,1)*ones(2*M,1);
    cutContour = ones(2*M,1);
end

% Create an array that gives the cut index for each vertex
vrt2cut = (1:M) .* [1; 1];
vrt2cut = vrt2cut(:);

% Determine what vertices are adjacenct to each other
pdistInds = getPdistInds(2*M);

% Get absolute indices difference between pairs
d = abs(diff(cutInds(pdistInds),[],2));

% Determine pairs from same contour
validPairs = cutContour(pdistInds(:,1)) == cutContour(pdistInds(:,2));

% Get contour length of each pair
N = contourLength(pdistInds(:,1));

% Determine if the distance is closer directly or wrapped around
upDown = round(d ./ N);

% Wrap the distance if necessary
d = (upDown .* N) + ((-1).^upDown) .* d;

% Get the vertices that are adjacent
% [d, validPairs]

if numel(ADJACENCY_DISTANCE) > 1
    cutA = cell(1,numel(ADJACENCY_DISTANCE));
    vrtA = cell(1,numel(ADJACENCY_DISTANCE));
    
    for i = 1:numel(ADJACENCY_DISTANCE)
        adjacentVrts = pdistInds((d <= ADJACENCY_DISTANCE(i)) & validPairs,:);

        % Construct the (symmetric) adjacency matrix

        adjacentCuts = vrt2cut(adjacentVrts);
        if size(adjacentVrts,1)==1
            adjacentCuts = adjacentCuts';
        end

        vrtA{i} = accumarray(adjacentVrts,1,[2*M,2*M]);
        vrtA{i} = vrtA{i} + vrtA{i}';

        cutA{i} = accumarray(adjacentCuts,1,[M,M]);
        cutA{i} = cutA{i} + cutA{i}';
    end
    
else
    adjacentVrts = pdistInds((d <= ADJACENCY_DISTANCE) & validPairs,:);

    % Construct the (symmetric) adjacency matrix

    adjacentCuts = vrt2cut(adjacentVrts);
    if size(adjacentVrts,1)==1
        adjacentCuts = adjacentCuts';
    end

    vrtA = accumarray(adjacentVrts,1,[2*M,2*M]);
    vrtA = vrtA + vrtA';

    cutA = accumarray(adjacentCuts,1,[M,M]);
    cutA = cutA + cutA';
end

end
% computeCutAdjacency : changeLog
% 2016-10-21 : added vertex adjacency to output.
% 2016-10-23 : fixed but with nan cutInds (which are from triangle centers)
% 2016-10-27 : allow adjacencyDistance to have more than one value.
% switched the comparison with adjacnecyDistance to <= instead of <



% function A = computeCutAdjacency__vertex(cutsVrts)
% % COMPUTECUTADJACENCY  Determine which cuts are adjacent to each other and
% % return the adjacency matrix.
% %
% % Two cuts are adjacent if they have vertices within ADJACENCY_DISTANCE of
% % each other
% 
% % James Kapaldo
% % 2016-10-20
% 
% ADJACENCY_DISTANCE = 3;
% 
% M = size(cutsVrts,1);
% 
% % Create an array that gives the cut index for each vertex
% vrt2cut = (1:M) .* [1; 1];
% vrt2cut = vrt2cut(:);
% 
% % Get the distance between the cut vertices
% cutsVrts = reshape(cutsVrts',2,2*M)';
% 
% % Determine what vertices are adjacenct to each other
% pdistInds = getPdistInds(2*M);
% d = pdistmex(cutsVrts','euc',[])';
% adjacentVrts = pdistInds(d <= ADJACENCY_DISTANCE,:);
% 
% % Construct the (symmetric) adjacency matrix
% A = accumarray(vrt2cut(adjacentVrts),1,[M,M]);
% A = A + A';
% 
% end