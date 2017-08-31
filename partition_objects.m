function [BW_partitioned, cuts, Info] = partition_objects(I,BW,r0,options)
% PARTITION_OBJECTS Declump the nuclei in an image
%
% [BW,cuts,Info] = computeNucleiCenters(I,BW,options)
%
% I - input image
% BW - object mask for image
% r0 - cell array of seed-points for each object
% options - declumpOptions class object
%
% BW - object mask after partitioning clumps
% cuts - Nx4 array. cuts(i,1:2) and cuts(i,3:4) give the two vertices of
%        the i'th cut
% Info - cell array of structures giving information

% James Kapaldo



% Get number of rows in image
nRows = size(BW,1);

% Remove small holes from mask
BW = bwareaopen(BW, options.Minimum_Hole_Size);
BW = ~bwareaopen(~BW, options.Minimum_Hole_Size,4);
CC = bwconncomp(BW); % Get connected commponents
pixelList = CC.PixelIdxList;

% Get object scales
objectScale = compute_objectScale(BW, pixelList);

% Create image and image edges arrays -----------------------------------
I = imfilter(I, fspecial('gaussian',7,1)); % Smooth the image
S = getImageEdges(I,options.Use_GPU);
I = cellfun(@(x) mean(I(x))./I(x), pixelList,'UniformOutput',false); % Slice 1 over image intensity
S = cellfun(@(x) S(x), pixelList,'UniformOutput',false); % Slice image edges

% Compute boundary information
[B,N,K] = computeBoundaryInformation(BW, objectScale, options);

% If there is an object of interest, remove all the others.
[pixelList, B, N, K, I, S, r0, objectScale] = getObjectOfInterest(options.Object_Of_Interest, pixelList, B, N, K, I, S, r0, objectScale);

% Do not partition if the object is convex
isConvex = cellfun(@(x,y) ~any(x(~isnan(x)) > (0.25./y)),K, num2cell(objectScale)); % object convex?

% Offset the r0set points to coorespond to object origin
objOffset = cellfun(@(x) min(x,[],1,'omitnan') - 1, B,'UniformOutput',false);
r0 = cellfun(@(x,y) x - y, r0, objOffset,'UniformOutput',false);
B = cellfun(@(x,y) x - y, B, objOffset,'UniformOutput',false);

% Compute the seed points.
% Note, if the seed-points also need to be calculated then it could be good
% to put that calculation (also add any needed inputs).
[cuts, Info] = processObjects(pixelList, B, N, K, I, S, r0, isConvex, nRows, options);

toRemove = cellfun(@isempty, cuts);
objOffset(toRemove) = [];
cuts(toRemove) = [];

if isempty(cuts)
    BW_partitioned = BW;
else

    % Offset seed points to image coordinates
    cuts = cellfun(@(x,y) x + [y, y], cuts, objOffset,'UniformOutput',false);

    % Offset r0 and r_final to image coordinates
    if options.Debug && ~isempty(objOffset)
        for i = 1:numel(Info)
            Info{i}.objOffset = objOffset{i};
        end
    end

    % Convert seedPoints from cell array to array, Nx3
    cuts = cat(1,cuts{:});

    BW_cuts = accumLines(cuts(:,1:2),cuts(:,3:4),size(BW),8);
    BW_partitioned = BW & ~BW_cuts;
    BW_partitioned = bwareaopen(BW_partitioned,options.Minimum_Hole_Size,4);
end

if options.Use_GPU
    gpuDevice([]);
end

if all(cellfun(@isempty, Info))
    Info = [];
end

end
