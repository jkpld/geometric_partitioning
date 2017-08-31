function [BW,cuts,Info,intraSS] = declumpNuclei(I,BW,options)
% DECLUMPNUCLEI Declump the nuclei in an image
%
% [BW,cuts,Info] = declumpNuclei(I,BW,options)
%
% I - input image
% BW - object mask for image
% options - declumpOptions class object
%
% BW - object mask after partitioning clumps
% cuts - Nx4 array. cuts(i,1:2) and cuts(i,3:4) give the two vertices of
%        the i'th cut
% Info - cell array of structures giving information

% James Kapaldo

% if options.Use_GPU
%     dev = gpuDevice();
% end

% Get number of rows in image
numImRows = size(I,1);

% Smooth the image
Is = imfilter(I,fspecial('gaussian',7,1));

% Get image edges
S = getImageEdges(Is,options.Use_GPU);

% Remove small holes from mask
BW = ~bwareaopen(~BW,options.Minimum_Hole_Size,4);

% Get connected commponents
CC = bwconncomp(BW);

% Get intra object edges
intraSS = intraObjectEdges(Is,BW,CC,options);

% Remove small holes from mask
BW = ~bwareaopen(~BW,options.Minimum_Hole_Size,4);

% Get connected commponents
CC = bwconncomp(BW);
pixelList = CC.PixelIdxList;

% Create sliced cell arrays with image edges, and 1/image intensities
Is = double(Is);
S = cellfun(@(x) S(x), pixelList,'UniformOutput',false);
Iunder1 = cellfun(@(x) mean(Is(x))./Is(x), pixelList,'UniformOutput',false);
intraS = cellfun(@(x) intraSS(x), pixelList,'UniformOutput',false);

% Compute boundary information
[B,N,K,M] = computeBoundaryInformation(BW,options);


% If there is an object of interest, remove all the others.
if ~isempty(options.Object_Of_Interest)
    B = B(options.Object_Of_Interest);
    N = N(options.Object_Of_Interest);
    K = K(options.Object_Of_Interest);
    M = M(options.Object_Of_Interest);
    pixelList = pixelList(options.Object_Of_Interest);
    S = S(options.Object_Of_Interest);
    Iunder1 = Iunder1(options.Object_Of_Interest);
    intraS = intraS(options.Object_Of_Interest);
    options.Use_Parallel = false;
end

% Initialize sliced variables
cutPixList = cell(1,numel(B));
cuts = cell(1,numel(B));
Info = cell(1,numel(B));

% Determine if debugging plots should be produced
produceDebugPlot = false;
if numel(B) == 1 && options.Debug && ~options.Use_Parallel
    produceDebugPlot = true;
end

% Convert the options to a normal structure so that we dont have to copy it
% and re-intialize each copy.
warning('off','MATLAB:structOnObject')
options = struct(options);
warning('on','MATLAB:structOnObject')

% Declump objects

if options.Use_Parallel
    
    parfor obj = 1:numel(B)
    
        [cutPixList{obj},cuts{obj},Info{obj}] = declumpObject__(B{obj},N{obj},K{obj},M{obj},S{obj},Iunder1{obj},intraS{obj},pixelList{obj},numImRows,options,produceDebugPlot);

        if ~isempty(Info{obj}.error)
            % If there was an error, try again as most error occure
            % because the centers were in just the wrong position and
            % they will not be again.
            [cutPixList{obj},cuts{obj},Info{obj}] = declumpObject__(B{obj},N{obj},K{obj},M{obj},S{obj},Iunder1{obj},intraS{obj},pixelList{obj},numImRows,options,produceDebugPlot);

            if ~isempty(Info{obj}.error)
                fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', obj, obj)
                fprintf(2,'%s\n', getReport(Info{obj}.error,'basic'))
                fprintf('\n')
            end % if
        end % if
    end % parfor
else
    
    generateDisplayAt = unique(round(linspace(1,numel(B),7)));
    processTimes = zeros(1,numel(B));
    
    for obj = 1:numel(B)
        procTime = tic;
        
        [cutPixList{obj},cuts{obj},Info{obj}] = declumpObject__(B{obj},N{obj},K{obj},M{obj},S{obj},Iunder1{obj},intraS{obj},pixelList{obj},numImRows,options,produceDebugPlot);

        if ~isempty(Info{obj}.error)
            % If there was an error, try again as most error occure
            % because the centers were in just the wrong position and
            % they will not be again.
            [cutPixList{obj},cuts{obj},Info{obj}] = declumpObject__(B{obj},N{obj},K{obj},M{obj},S{obj},Iunder1{obj},intraS{obj},pixelList{obj},numImRows,options,produceDebugPlot);

            if ~isempty(Info{obj}.error)
                fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', obj, obj)
                fprintf(2,'%s\n', getReport(Info{obj}.error,'basic'))
                fprintf('\n')
            end % if
        end % if
        
        processTimes(obj) = toc(procTime);
        if any(obj == generateDisplayAt)
            fprintf('%s >> %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),obj,numel(B),sum(processTimes(1:obj))/60,mean(processTimes(1:obj))*numel(B)/60)
        end % if
    end % for
end % if

if options.Use_GPU
    resetGPU
end

% Produce output
cutPixList = cat(1,cutPixList{:});
cuts = cat(1,cuts{:});

BW(cutPixList) = false;
BW = bwareaopen(BW,options.Minimum_Hole_Size);

if all(cellfun(@isempty, Info))
    Info = [];
end

end

function resetGPU
%     reset(dev);
    gpuDevice([]);
end



function [cutPixList, cuts, Info] = declumpObject__(B,N,K,M,S,I,intraS,pixList,numImRows,options,produceDebugPlot)


% Note
% Put the whole thing in a try-catch statement. As the particle positions
% can be in different places every single time, sometimes there are errors
% in the partitioning that have never been seen before (and thus not
% fixed). If there is an error in an object, then cutPixList and cuts will
% be empty. Info will have a field named 'error' that contains the error
% message identifier. Info will also have a field named
% 'centers_ObjectCoordinates' that contains the centers from the particle
% simulation in the objects coordinates (not the image coordinates). To
% debug the error, hard code in these center values into current_centers
% below the call to computeObjectCenters() below.


try


% Initialize output
cutPixList = [];
cuts = [];
if options.Debug
    Info = struct('r0',[],'r_end',[],'solverTime',NaN,'centers',[],'convexOrTooSmall',true,'cutCalculationTime',NaN,'error',[]);
else
    Info.error = [];
end
current_centers = [];

% Create three small images with the mask, the edges, and 1 over the image
% intensity.
[BW,S,I,intraS] = createObjectImages(pixList,S,I,intraS,numImRows);

% Area of object
objArea = numel(pixList);

% Only run declumping if there are concave regions and the region is larger
% than our minimum size.
if ~any(K>0) && (objArea < (pi*options.Wigner_Seitz_Radius^2))
    return;
end

% Compute the top left position of the boundary and offset the
% boundary and markers.
topLeftB = min(B,[],1,'omitnan') - 1;
B = B - topLeftB;
M = M - topLeftB;

% Compute the centers of the object
if options.Debug
    
    [current_centers,Info] = computeObjectCenters(BW,B,intraS,M,options);
    
    if isempty(current_centers)
        Info.r_end = [];
        Info.centers = [];
    else
        Info.r_end = Info.r_end + topLeftB;
        Info.centers = current_centers + topLeftB;
    end
    Info.r0 = Info.r0 + topLeftB;
    Info.topLeft = topLeftB;
    Info.convexOrTooSmall = false;
    Info.cutCalculationTime = NaN;
    Info.error = [];
else
    
    current_centers = computeObjectCenters(BW,B,intraS,M,options);
    
end

% current_centers = ...
%    [16.0000   24.2174;
%    71.8846   33.8846];
% 
% figure(1)
% clf(1)
% line(B(:,1),B(:,2))
% hold on
% plot(current_centers(:,1),current_centers(:,2),'ro')

% If there is not more than one center, then duclumping will not be
% run on this object, continue on to next.
if size(current_centers,1)<2
    Info.cutCalculationTime = NaN;
    return;
end

% Start a clock for timing.
if options.Debug
    start1 = tic;
end

% Compute the object paritions
if options.Debug
    [objCuts,~,cutCenters,cutIndices,CreateObjectPartitionInfo] = createObjectPartitions(B,N,current_centers,options);
    Info.CreateObjectPartition = CreateObjectPartitionInfo;
    Info.originalCuts = objCuts + [topLeftB, topLeftB];
    
    if produceDebugPlot
        [fig,ax] = plotObjectCutDebugInfo(B,N,objCuts,Info);
        ax.Title.String = 'Cuts before optimization';
        drawnow;
        goDark(fig)
    end
else
    [objCuts,~,cutCenters,cutIndices] = createObjectPartitions(B,N,current_centers,options);
end

% Create the triangle groups
triangleGroups = constrainedObjectCenterTriangulation(B,current_centers,options);

% Get cut association before optimization
[~,vrtA] = computeCutAdjacency(cutIndices,B,5);
associatedCuts = findAssociatedCuts(triangleGroups,objCuts,current_centers,vrtA,cutCenters);

% Optimize the cuts
[objCuts,objCutNs,objCutKs,cutIndices] = optimizeCuts(cutIndices,B,N,K,options);

if options.Debug
    Info.optimizedCuts = objCuts + [topLeftB, topLeftB];
    Info.cutsIntersectBoundary = checkCutIntegrity(objCuts,B,options.Use_GPU);
    
    if produceDebugPlot
        [fig,ax] = plotObjectCutDebugInfo(B,N,objCuts,Info);
        ax.Title.String = 'Cuts after optimization';
        drawnow;
        goDark(fig)
    end
end

% Get cut association after optimization
[~,vrtA] = computeCutAdjacency(cutIndices,B,5);
associatedCutsAfterOptimization = findAssociatedCuts(triangleGroups,objCuts,current_centers,vrtA,cutCenters);

% Combine associations before and after to get the complete
% list.
for i = 1:numel(associatedCuts)
    associatedCuts{i} = unique([associatedCuts{i};associatedCutsAfterOptimization{i}]);
end

% Create triangle cuts
[objTriCuts,objTriCutNs,objTriCutKs,~,InfoTri] = createTriangleCuts(B, N, K, triangleGroups, current_centers, associatedCuts, options);

if options.Debug
    
    if ~isempty(objTriCuts)
        for i=1:numel(InfoTri)
            InfoTri(i).originalCuts = InfoTri(i).originalCuts + [topLeftB, topLeftB];
            InfoTri(i).vertexOptimizedCuts = InfoTri(i).vertexOptimizedCuts + [topLeftB, topLeftB];
            InfoTri(i).optimizedCuts = InfoTri(i).optimizedCuts + [topLeftB, topLeftB];
            InfoTri(i).triangleCenterSearchPoints = InfoTri(i).triangleCenterSearchPoints + topLeftB;
        end
    end
    Info.CreateTriangleCuts = InfoTri;
    if produceDebugPlot
        if isempty(triangleGroups)
            fprintf('No triangle information to plot!\n')
        else
            [fig,ax] = plotObjectTriDebugInfo(B,objTriCuts,objTriCutKs,Info);
            ax.Title.String = 'Triangle cuts';
            drawnow;
            goDark(fig)
        end
    end
end

% Determine winning cuts
[~,vrtA] = computeCutAdjacency(cutIndices,B,1);
[objCuts,~,objCutKs,isTriCut,scores] = chooseWinningCuts(objCuts,objCutNs,objCutKs,objTriCuts,objTriCutNs,objTriCutKs,associatedCuts,vrtA,S,I);

if options.Debug
    Info.cutScores = scores;
end

% Remove any cuts that are not at concave points.
if ~isempty(objCuts) && any(~isTriCut)
%     tmpK = objCutKs;
%     tmpK(isnan(tmpK)) = inf;
%     toRemove = ~all(tmpK < 1/options.Max_Radius,2);
    toRemove = (mean(objCutKs,2) < 1/(2*options.Max_Radius)) & ~isTriCut; % maybe use max instead of mean.
    objCuts(toRemove,:) = [];
end

% Done!

% Offset the cuts by the topLeftB
objCuts(:,1:2) = objCuts(:,1:2) + topLeftB;
objCuts(:,3:4) = objCuts(:,3:4) + topLeftB;

% Get the linear indices of the pixels for each cut.
cutPixList_cell = cell(size(objCuts,1),1);
for current_cut = 1:size(objCuts,1)
    inds = rayTrace(objCuts(current_cut,1:2),objCuts(current_cut,3:4));
    cutPixList_cell{current_cut} = inds(:,1) + (inds(:,2)-1)*numImRows;
end

cutPixList = cat(1,cutPixList_cell{:});
cuts = objCuts;
if options.Debug
    Info.cutCalculationTime = toc(start1);
end


catch ME
   cutPixList = [];
   cuts = [];
   Info.error = ME;
   Info.centers_ObjectCoordinates = current_centers;
end


end



function [fig,ax] = plotObjectCutDebugInfo(B,N,cuts,Info)
% PLOTOBJECTCUTDEBUGINFO  Helper function for plotting debuggin information

% James Kapaldo
% 2016-10-29

centers = Info.centers - Info.topLeft;
r0 = Info.r0 - Info.topLeft;
r_end = Info.r_end - Info.topLeft;

fig = figure;
ax = axes('Parent',fig);
hold on

aC = Info.CreateObjectPartition.assignedCenter;
assignmentLine = nan(3*size(B,1),2);
ind = 1:3:3*size(B,1);
assignmentLine(ind,:) = B;
assignmentLine(ind(~isnan(aC))+1,:) = centers(aC(~isnan(aC)),:);

p1 = plot(assignmentLine(:,1), assignmentLine(:,2),'Color',0.4*[1 1 1]);
plot(B(:,1), B(:,2),'Color','k','LineWidth',2)
p2 = plot(centers(:,1),centers(:,2),'bo','MarkerFaceColor','b');
p3 = quiver(B(:,1),B(:,2),N(:,1),N(:,2),0.25,'Color',0.5*[1 1 1],'LineWidth',1.5);
p4 = plot(r0(:,1),r0(:,2),'r.');
p5 = plot(r_end(:,1),r_end(:,2),'b.');

for i = 1:numel(Info.CreateObjectPartition.CalculateCuts.partitions)
    pP = line(Info.CreateObjectPartition.CalculateCuts.partitions{i}(:,1),Info.CreateObjectPartition.CalculateCuts.partitions{i}(:,2),'color','g');
end


for i = 1:size(cuts,1)
    tmp = (cuts(i,1:2) + cuts(i,3:4))/2;
    pC = line(cuts(i,[1,3]),cuts(i,[2,4]),'color',[0 0 1],'lineWidth',2);
    text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color',[0 0 0.5],'Tag','ignore','Clipping','on')
end

p6 = plot(cuts(:,1),cuts(:,2),'o','LineWidth',1,'Color',0.2*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g');
p7 = plot(cuts(:,3),cuts(:,4),'d','LineWidth',1,'Color',0.2*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g');

legend([p1,p3,p4,p5,p2,pP,pC,p6,p7],{'Boundary assignments','Boundary normals','Initial particle locations','Final particle locations','Centers','Partitions before optimization','Cuts','Cut starts','Cut ends'},'Location','eastoutside')

for i = 1:size(centers,1)
    text(centers(i,1),centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color',[0.5 0 0],'Tag','ignore','Clipping','on')
end

grid on
box on
axis tight
daspect([1 1 1])
drawnow;


end


function [fig,ax] = plotObjectTriDebugInfo(B,cut,cutK,Info)
% PLOTOBJECTTRIDEBUGINFO  Helper function for plotting debuggin information
% about triangle cuts

% James Kapaldo
% 2016-10-29

centers = Info.centers - Info.topLeft;
TriInfo = Info.CreateTriangleCuts;

fig = figure;
ax = axes('Parent',fig);
hold on

plot(B(:,1), B(:,2),'Color','k','LineWidth',2)
pC = plot(centers(:,1),centers(:,2),'bo','MarkerFaceColor','b');

for group = 1:numel(TriInfo)
    
    for i = 1:size(TriInfo(group).triGroup,2)
        triCent = TriInfo(group).triGroup(:,i);
        p4 = plot(centers(triCent([1,2]),1),centers(triCent([1,2]),2),'Color',0.3*[1 1 1],'LineWidth',2);
        plot(centers(triCent([1,3]),1),centers(triCent([1,3]),2),'Color',0.3*[1 1 1],'LineWidth',2);
        plot(centers(triCent([2,3]),1),centers(triCent([2,3]),2),'Color',0.3*[1 1 1],'LineWidth',2);
    end
    
    pM = plot(TriInfo(group).triangleCenterSearchPoints(:,1)-Info.topLeft(:,1),TriInfo(group).triangleCenterSearchPoints(:,2)-Info.topLeft(:,2),'ys','MarkerSize',4,'LineWidth',0.5);
    for i = 1:size(TriInfo(group).originalCuts)
        p1 = plot(TriInfo(group).originalCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).originalCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[1 0 0],'LineWidth',2);
        p2 = plot(TriInfo(group).vertexOptimizedCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).vertexOptimizedCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[0 1 0],'LineWidth',2);
        p3 = plot(TriInfo(group).optimizedCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).optimizedCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[0 0 1],'LineWidth',2);
    end
    
    cutG = cut{group}(isnan(cutK{group}(:,1)),:);
    for i = 1:size(cutG,1)
        plot(cutG(i,[1,3]),cutG(i,[2,4]),'Color',0.5*[0 0 1],'LineWidth',2,'LineStyle','--')
    end
end

legend([pC,p4,pM,p1,p2,p3],{'Centers','Triangle edges','Triangle center search points','Original cuts','Vertex optimized cuts','OptimizedCuts'},'Location','eastoutside')

for i = 1:size(centers,1)
    text(centers(i,1),centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color',[0.5 0 0],'Tag','ignore','Clipping','on')
end

grid on
box on
axis tight
daspect([1 1 1])
drawnow;

end

