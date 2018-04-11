%% Example Usage with SALR Clustering
% Use geometric partitioning to partition clumps of nuclei using SALR
% clustering to find the seed-points.
%
% Note, must have SALR_Clustering installed.

%% Load image and mask
I = imread('image_67.tif');
BW = logical(imread('mask_67.tif'));

% Set partitioning options and SALR options
partOptions = partitionOptions();
partOptions.Use_GPU = 1;
partOptions.Minimum_Hole_Size = 10;

salrOptions = seedPointOptions();
salrOptions.Wigner_Seitz_Radius = 5;
salrOptions.Point_Selection_Method = 'r0set_uniformRandom';
salrOptions.Maximum_Initial_Potential = 1/5;
salrOptions.Max_Distance_Transform = 18;
salrOptions.Potential_Parameters = [-1, 2, 13];
salrOptions.Verbose = true;
salrOptions.Use_Parallel = true;

% Find seed-points and partition nuclei clumps
[BW_part, cuts, seedPoints, Info] = declumpNuclei(I, bwconncomp(BW), salrOptions, partOptions);

% Remove small objects
BW_part = bwareaopen(BW_part,partOptions.Minimum_Object_Size);

%% Plot results
plot_partition_results(BW_part, seedPoints, [700, 1970; 1450, 2670])