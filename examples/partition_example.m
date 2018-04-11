%% Example Usage
% Use geometric partitioning to partition clumps of nuclei using the true
% nuclei centers.

% Load image, mask, and seed-points
I = imread('image_67.tif');
BW = imread('mask_67.tif');
r0 = load_truth('67');

% Create the partition options

partOptions = partitionOptions();
partOptions.Minimum_Hole_Size = 10;
partOptions.Use_GPU = 1;

% Partition objects 
[BW_partitioned, cuts, Info] = partition_objects(I, BW, r0, partOptions);

% Plot results
plot_partition_results(BW_partitioned, r0, [700, 1970; 1450, 2670])
