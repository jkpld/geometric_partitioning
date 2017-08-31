function [cuts, Info] = geometric_partition(B, N, K, I, S, r0, options)
% GEOMETRIC_PARTITION
%
% Input parameters -
%   B : n x 2 array of boundary vertices
%   N : n x 2 array of boundary vertex normals
%   K : n x 1 array of boundary curvature
%   I : image of object
%   S : edges of I
%   r0 : m x 2 array of seed points
%   options : ...needed options
%
% Output parameters -
%   cuts : w x 4 array of the end points of each cut
%
% B and r0 should be in the same coordinate system where the top left
% corner of the boundary is [1, 1]

cuts = [];
Info = [];

DEBUG = options.Debug;
DEBUG_PLOT = options.Debug_plot;


% Compute the object paritions
if DEBUG
    timer_start = tic;
    [objCuts,~,cutCenters,cutIndices,CreateObjectPartitionInfo] = createObjectPartitions(B,N,r0,options);
    Info.CreateObjectPartition = CreateObjectPartitionInfo;
    Info.originalCuts = objCuts;
    Info.centers = r0;
    if DEBUG_PLOT
        [fig,ax] = plotObjectCutDebugInfo(B,N,objCuts,Info);
        ax.Title.String = 'Cuts before optimization';
        drawnow;
        setTheme(fig,'dark')
    end
else
    [objCuts,~,cutCenters,cutIndices] = createObjectPartitions(B,N,r0,options);
end

if isempty(objCuts)
    return;
end

% Create the triangle groups
triangleGroups = constrainedObjectCenterTriangulation(B,r0,options);

% Get cut association before optimization
[~,vrtA] = computeCutAdjacency(cutIndices,B,5);
associatedCuts = findAssociatedCuts(triangleGroups,objCuts,r0,vrtA,cutCenters);

% Optimize the cuts
[objCuts,objCutNs,objCutKs,cutIndices] = optimizeCuts(cutIndices,B,N,K,options);

if DEBUG
    Info.optimizedCuts = objCuts;
    Info.cutsIntersectBoundary = checkCutIntegrity(objCuts,B,options.Use_GPU);
    
    if DEBUG_PLOT
        [fig,ax] = plotObjectCutDebugInfo(B,N,objCuts,Info);
        ax.Title.String = 'Cuts after optimization';
        drawnow;
        setTheme(fig,'dark')
    end
end

% Get cut association after optimization
[~,vrtA] = computeCutAdjacency(cutIndices,B,5);
associatedCutsAfterOptimization = findAssociatedCuts(triangleGroups,objCuts,r0,vrtA,cutCenters);

% Combine associations before and after to get the complete
% list.
for i = 1:numel(associatedCuts)
    associatedCuts{i} = unique([associatedCuts{i};associatedCutsAfterOptimization{i}]);
end

% Create triangle cuts
[objTriCuts,objTriCutNs,objTriCutKs,~,InfoTri] = createTriangleCuts(B, N, K, triangleGroups, r0, associatedCuts, options);

if DEBUG
    
%     if ~isempty(objTriCuts)
%         for i=1:numel(InfoTri)
%             InfoTri(i).originalCuts = InfoTri(i).originalCuts + [topLeftB, topLeftB];
%             InfoTri(i).vertexOptimizedCuts = InfoTri(i).vertexOptimizedCuts + [topLeftB, topLeftB];
%             InfoTri(i).optimizedCuts = InfoTri(i).optimizedCuts + [topLeftB, topLeftB];
%             InfoTri(i).triangleCenterSearchPoints = InfoTri(i).triangleCenterSearchPoints + topLeftB;
%         end
%     end
    Info.CreateTriangleCuts = InfoTri;
    if DEBUG_PLOT
        if isempty(triangleGroups)
            fprintf('No triangle information to plot!\n')
        else
            [fig,ax] = plotObjectTriDebugInfo(B,objTriCuts,objTriCutKs,Info);
            ax.Title.String = 'Triangle cuts';
            drawnow;
            setTheme(fig,'dark');
        end
    end
end

% Determine winning cuts
[~,vrtA] = computeCutAdjacency(cutIndices,B,1);
[objCuts,~,objCutKs,isTriCut,scores] = chooseWinningCuts(objCuts,objCutNs,objCutKs,objTriCuts,objTriCutNs,objTriCutKs,associatedCuts,vrtA,S,I);

if DEBUG
    Info.cutScores = scores;
end

% Remove any cuts that are not at concave points.
if ~isempty(objCuts) && any(~isTriCut)
%     tmpK = objCutKs;
%     tmpK(isnan(tmpK)) = inf;
%     toRemove = ~all(tmpK < 1/options.Max_Radius,2);
    toRemove = (mean(objCutKs,2) < 1/(2*options.Max_Radius)) & ~isTriCut; % maybe use max instead of mean.
%     toRemove = any(objCutKs < 1/(2*options.Max_Radius),2) & ~isTriCut;
    objCuts(toRemove,:) = [];
end

% Done!

% Offset the cuts by the topLeftB
% objCuts(:,1:2) = objCuts(:,1:2) + topLeftB;
% objCuts(:,3:4) = objCuts(:,3:4) + topLeftB;
% 
% % Get the linear indices of the pixels for each cut.
% cutPixList_cell = cell(size(objCuts,1),1);
% for current_cut = 1:size(objCuts,1)
%     inds = rayTrace(objCuts(current_cut,1:2),objCuts(current_cut,3:4));
%     cutPixList_cell{current_cut} = inds(:,1) + (inds(:,2)-1)*numImRows;
% end
% 
% cutPixList = cat(1,cutPixList_cell{:});
cuts = objCuts;

if DEBUG
    Info.cutCalculationTime = toc(timer_start);
end

