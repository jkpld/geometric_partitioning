function [BW_partitioned, cuts, seedPoints, Info] = declumpNuclei(I,CC,salr_options,geoPart_options,attempt)
% DECLUMPNUCLEI Use SALR clustering to find the nuclei centers in a binary
% image of nuclei and then partition the nuclei clumps using geometric
% partitioning.
%
% [seedPoints, Info] = declumpNuclei(I, BW, salr_options, geoPart_options)
% [seedPoints, Info] = declumpNuclei(I, BW, salr_options, geoPart_options, 'Object_Of_Interest', OoI)
%
% For each object in the image, compute the centers of curvature to use as
% initial points, then find the seed-points for each object.
%
% Input parameters:
% I : Nuclei image.
% CC : connected compontents structure for each object..
% salr_options : An instance of class seedPointOptions.
% geoPart_options : An instance of class partitionOptions.
% attempt : (optional) Nx1 logical array, where N is the number of objects,
%   that determines if partitioning should be attempted.
%
% Output parameters:
% BW_partitioned : Binary mask after partioning the nuclei
% cuts : Nx4 array. cuts(i,1:2) and cuts(i,3:4) give the two vertices of
%   the i'th cut.
% seedPoints : Nx3 array where the first two columns give the x,y location
%   of a nuclei center and the third column gives the object number to
%   which the nuclei belongs.
% Info : Cell array with a length equal to the number of objects. Each
%   element is a structure with 3 fields
%       The Info structure returned by computeObjectSeedPoints().
%       The Info structure returned by geometric_partition().
%       The object offset. (The upper left corner of the object in the
%         image.)

% There is no error checking on the inputs!

if nargin < 5
    attempt = true(CC.NumObjects,1);
end

% Extract pixels
pixelList = CC.PixelIdxList;

% Get number of rows in image
nRows = size(I,1);

% Form binary mask
BW = false(size(I));
BW(cat(1,pixelList{:})) = true;

% Get object scales
objectScale = compute_objectScale(BW, pixelList);

% Create image and image edge arrays -----------------------------------
I = imfilter(I, fspecial('gaussian',7,1)); % Smooth the image
S = getImageEdges(I,geoPart_options.Use_GPU);
I = cellfun(@(x) mean(I(x))./I(x), pixelList,'UniformOutput',false); % Slice 1 over image intensity
S = cellfun(@(x) S(x), pixelList,'UniformOutput',false); % Slice image edges

% Compute boundary information
[B,N,K,r0set] = computeBoundaryInformation(BW, objectScale, salr_options);

% If there is an object of interest, remove all the others.
% [pixelList, B, N, K, I, S, r0set, objectScale] = getObjectOfInterest(Object_Of_Interest, pixelList, B, N, K, I, S, r0set, objectScale);

% Use the object centroid as the seed point for any object that is convex
% or smaller than the particle area
area = cellfun(@numel,pixelList); % object areas
isConvex = cellfun(@(x,y) ~any(x(~isnan(x)) > (0.25./y)),K, num2cell(objectScale)); % object convex?
useCentroid = isConvex(:) | (area(:) < pi*salr_options.Wigner_Seitz_Radius.^2) | ~attempt(:);


% Offset the r0set points to coorespond to object origin
objOffset = cellfun(@(x) min(x,[],1,'omitnan') - 1, B,'UniformOutput',false);
r0set = cellfun(@(x,y) x - y, r0set, objOffset,'UniformOutput',false);
B = cellfun(@(x,y) x - y, B, objOffset,'UniformOutput',false);

% Compute the seed points.
n = numel(pixelList);
verbose = salr_options.Verbose;
Use_Parallel = salr_options.Use_Parallel;
if n == 1
    Use_Parallel = 0;
end

% Initialize sliced variables
Info = cell(n,1);
seedPoints = cell(n,1);
cuts = cell(n,1);

% Create a progress monitor
progres = displayProgress(n, 10, verbose, Use_Parallel, 'name', 'Partitioning objects, ');
Que = progres.start();

if Use_Parallel
    % If we are computing in parallel, then first convert the options class
    % element to a structure to prevent reinitiallization on transfer to
    % each worker.
    warning('off','MATLAB:structOnObject')
    salr_options = struct(salr_options);
    geoPart_options = struct(geoPart_options);
    warning('on','MATLAB:structOnObject')

    parfor obj = 1:n
        [cuts{obj}, seedPoints{obj}, Info{obj}] = processObject(pixelList{obj}, nRows, r0set{obj}, useCentroid(obj), B{obj}, N{obj}, K{obj}, I{obj}, S{obj}, obj, salr_options, geoPart_options);
        if ~isempty(Que), send(Que, obj), end
    end
else
    for obj = 1:n
        [cuts{obj}, seedPoints{obj}, Info{obj}] = processObject(pixelList{obj}, nRows, r0set{obj}, useCentroid(obj), B{obj}, N{obj}, K{obj}, I{obj}, S{obj}, obj, salr_options, geoPart_options);
        progres.iteration_end()
    end
end

delete(progres)

% Offset seed points to image coordinates
seedPoints = cellfun(@(x,y) [x(:,1:2) + y, x(:,3)], seedPoints, objOffset,'UniformOutput',false);

% Convert seedPoints from cell array to array, Nx3
seedPoints = cat(1,seedPoints{:});

% Save object offsets if debugging
if salr_options.Debug
    for i = 1:numel(Info)
        Info{i}.objOffset = objOffset{i};
    end
end

% Remove any empty cuts
toRemove = cellfun(@isempty, cuts);
objOffset(toRemove) = [];
cuts(toRemove) = [];

if isempty(cuts)
    BW_partitioned = BW;
else
    % Offset cuts to image coordinates
    cuts = cellfun(@(x,y) x + [y, y], cuts, objOffset,'UniformOutput',false);

    % Convert cuts from cell array to array, Nx4
    cuts = cat(1,cuts{:});

    BW_cuts = accumLines(cuts(:,1:2),cuts(:,3:4),size(BW),8);
    BW_partitioned = BW & ~BW_cuts;
end



if all(cellfun(@isempty, Info))
    Info = [];
end

end

function [cuts, seedPoints, Info] = processObject(pixList, nRows, r0set, useCentroid, B, N, K, I, S, obj, salr_options, geoPart_options)


    % Create mask and interior potential images for the object
    [objBW, objI, objS] = createObjectImages(pixList, nRows, true(numel(pixList),1), I, S);
    
    % Compute seed points    
    [seedPoints, salr_Info] = computeObjectSeedPoints(logical(objBW), salr_options, 'r0set', r0set, 'useCentroid', useCentroid, 'objNumber', obj);
    
    % Partition clump
    if size(seedPoints,1) > 1
        [cuts, geoPart_Info] = geometric_partition(B, N, K, objI, objS, seedPoints, geoPart_options, obj);
    else
        cuts = [];
        geoPart_Info = [];
    end
    
    % Add object number as third column. (This is mostly just helpful when comparing against truth data, as the truth data is labeled by each object.)
    seedPoints = [seedPoints, obj*ones(size(seedPoints,1),1)];

    % Combine Info structures
    Info.salr_Info = salr_Info;
    Info.geoPart_Info = geoPart_Info;
    
end

%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
