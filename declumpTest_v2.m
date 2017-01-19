% function declumpTest_v2
% close all


% Parameters -------------------------------------------------------------

% Boundary Assignment
options.Max_Radius = 35;                                    % [pixels]
options.Min_Angle = 0.5;                                    % [dot product] : angle would be acos(Min_Angle)

% Object center calculation (particle simulation)
options.Wigner_Seitz_Radius = 10;                           % [pixels]
options.Initial_Speed = 0.01;                               % [pixels/time]
options.Point_Selection_Method = 'curvatureUniformRandom';  % {'Random','UniformRandom','curvatureRandom','curvatureUniformRandom'}
options.Potential_Depth = -1;                               % [arb. units.]
options.Potential_Minimum_Location = 2;                     % [pixels]
options.Potential_Extent = 15;                              % [pixels]

% Triangulation
options.Min_Interior_Angle = 20;                            % [degrees]
options.Max_Interior_Angle = 110;                           % [degrees]

% Optimization
options.Search_Radius = 7;                                  % [pixels]

% Mask
options.Minimum_Hole_Size = 50;                             % [pixels^2]

% Curvature
options.Curvature_Smoothing_Size = 2;                       % [pixels]

% Computation
options.Use_GPU = true;                                     % [boolean]
options.Use_Parallel = false;                               % [boolean]

% Debug : debuging outputs extra potentially useful information, this takes
% extra memory
options.Debug = true;                                       % [boolean]

%%
DEBUG = 1;

options.maxRadius = 35;
options.minAngle = 0.5;
options.searchRadius = 7;

% Works well for cancer
options.WIGNER_SEITZ_RADIUS = 10;
options.INITIAL_SPEED = 0.01;
options.pointSelectionMethod = 'curvatureUniformRandom';

depth = -1;
center = 2;
extent = 15;

% Works pretty well for OKF
% options.WIGNER_SEITZ_RADIUS = 6;
% options.INITIAL_SPEED = 0.01;
% options.pointSelectionMethod = 'curvatureUniformRandom';
% 
% 
% depth = -1;
% center = 2;
% extent = 9;

if ~exist('potentialParameters','var')
    potentialParameters = load('K:\Google_Drive\MATLAB\radiationLab\Projects\nucleiDeclumpByEM\potentialParameters.mat');
end

[~,i1] = min(abs(potentialParameters.depth-depth));
[~,i2] = min(abs(potentialParameters.center-center));
[~,i3] = min(abs(potentialParameters.extent-extent));

InteractionOptions.type = 'SRALRR'; % {'Coulomb', 'SRALRR'} % SRALRR = short-range attractive long-range repulsive
InteractionOptions.params = permute(potentialParameters.parameters(i1,i2,i3,:),[4,1,2,3]);

% params = squeeze(potentialParameters.parameters(i1,i2,:,:));

options.InteractionOptions = InteractionOptions;

KAPPA_SMOOTHING_SIGMA = 2;
MINIMUM_HOLE_AREA = 50;
% =======================================================================

rs = options.WIGNER_SEITZ_RADIUS;
minArea = pi*rs^2;
%
tic
% [I,Is,BW] = getDeclumpTestCase(6);
[I,Is,BW] = getDeclumpTestCase('6and7','cancer');
toc
%
% profile on
tic
S = getImageEdges(Is,false);
toc
% profile off
% profile viewer

% figure(10)
% clf(10)
% imshow(S,[])
% colorbar
% drawnow;
% goDark(gcf);

tic
BW = ~bwareaopen(~BW,MINIMUM_HOLE_AREA,4);
CC = bwconncomp(BW);
pixelList = CC.PixelIdxList;
S = cellfun(@(x) S(x), pixelList,'UniformOutput',false);
Iobjects = cellfun(@(x) I(x), pixelList,'UniformOutput',false);
toc
%%
% profile on
tic
[B,Bnorm,curvature,markers] = computeBoundaryInformation(BW,KAPPA_SMOOTHING_SIGMA,options.maxRadius,false); 
toc
% profile off
% profile viewer
%%
% dev = gpuDevice();

r0 = cell(numel(B),1);
r_end = cell(numel(B),1);
centers = cell(numel(B),1);
cuts = cell(numel(B),1);
Info = cell(numel(B),1);
areas = zeros(numel(B),1);
clc

sizeI1 = size(I,1);

offset_I = [0 0]; % Use this variable to store the offset of the top left pixel of image I from the very large image that I comes from.


start = tic;
% profile on
ind = 18*22+3;%9+15*22;
for obj = 1:numel(B)

    % Only run declumping if there are concave regions on the current
    % object.
    [objBW,objS,objI] = getMaskAndEdge(pixelList{obj},S{obj},Iobjects{obj},sizeI1);
%     objBW = createMaskFromBoundary(B{obj});
    areas(obj) = numel(pixelList{obj});
    
%     options.InteractionOptions.params = params(newExtentInd,:);
%     [~,i3] = min(abs(potentialParameters.extent-new_extent));
    
    if any(curvature{obj}>0) && areas(obj) > minArea
        try
            
        % Get the current object boundary, boundary' normals, and markers.    
        objB = B{obj};
        objN = Bnorm{obj};
        objK = curvature{obj};
        objMarkers = markers{obj};

        % Compute the top left position of the boundary and offset the
        % boundary and markers.
        topLeftB = min(objB,[],1,'omitnan') - 1;
        objB = objB - topLeftB;
        objMarkers = objMarkers - topLeftB;

        % *** Compute the centers of the object
        if DEBUG
            [current_centers,Info{obj}] = computeObjectCenters(objBW,objB,objMarkers,options);
        else
            current_centers = computeObjectCenters(objBW,objB,objMarkers,options);
        end

        if ~Info{obj}.converged
            fprintf('No converge in object %d\n',obj);
        end
%         current_centers(5,:) = [];
        % If there is not more than one center, then duclumping will not be
        % run on this object, continue on to next.
        if size(current_centers,1)>1

            start1 = tic;
            
            if DEBUG
                % Offset and save initial and final particle locations as
                % well as the centers.
                r0{obj} = Info{obj}.r0 + topLeftB;
                r_end{obj} = Info{obj}.r_end + topLeftB;
                centers{obj} = current_centers + topLeftB;
            end
        
            % *** Compute the object paritions
            if DEBUG
               [objCuts,~,cutCenters,cutIndices,CreateObjectPartitionInfo] = createObjectPartitions(B,N,current_centers,options);
                Info.CreateObjectPartition = CreateObjectPartitionInfo;
                
                if CreateObjectPartitionInfo.CalculateCuts.CleanCuts.unpairedCutRemoved
                    fprintf('duplicateCutRemoved in object %d\n',obj)
                end                
            else
                [objCuts,~,cutCenters,cutIndices] = createObjectPartitions(objB,objN,current_centers,options);
            end
            
            
            
%             objCuts
%             
%             figure(112)
%             clf(112)
% 
%             for i = 1:size(objB,1)
%                 if ~isnan(CreateObjectPartitionInfo.assignedCenter(i))
%                     line([objB(i,1),current_centers(CreateObjectPartitionInfo.assignedCenter(i),1)], [objB(i,2),current_centers(CreateObjectPartitionInfo.assignedCenter(i),2)],'color',0.6*[1 1 1])
%                 end
%             end
%             
%             line(objB(:,1),objB(:,2),'Color','k','LineWidth',2)
%             hold on
%             plot(current_centers(:,1),current_centers(:,2),'ro','MarkerFaceColor','r')
%             plot(objCuts(:,1),objCuts(:,2),'o','LineWidth',1,'Color',0.3*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g')
%             plot(objCuts(:,3),objCuts(:,4),'d','LineWidth',1,'Color',0.3*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g')
% %             
% %             line(r_end{obj}(:,1)-topLeftB(1),r_end{obj}(:,2)-topLeftB(2),'Marker','*','Color','y','LineStyle','none')
%             quiver(objB(:,1),objB(:,2),objN(:,1),objN(:,2),0.25,'Color',0.8*[1 1 1])
%             drawnow;
%             goDark(gcf)
%             for i = 1:size(objCuts,1)
%                 tmp = (objCuts(i,1:2) + objCuts(i,3:4))/2;
%                 line(objCuts(i,[1,3]),objCuts(i,[2,4]),'color',[0 0 1],'lineWidth',2)
%                 text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color','b')
%             end
%             
%             for i = 1:size(current_centers,1)
%                 text(current_centers(i,1),current_centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color','k')
%             end
%             
%             for i = 1:numel(CreateObjectPartitionInfo.partitions)
%                 line(CreateObjectPartitionInfo.partitions{i}(:,1),CreateObjectPartitionInfo.partitions{i}(:,2),'color','g')
%             end
% %             str = {'Note edge 4 should be connecting through the hole, but it is';
% %              'not because the normals at the hole and the boundary';
% %              'start/end points are perpendicular.'};
% %             title(str)
% %             grid on
% %             box on
% %             axis tight
%             daspect([1 1 1])
% %             error('somerror')
             
            
            % Create the triangle groups
            triangleGroups = constrainedObjectCenterTriangulation(objB,objN,current_centers);

            % Get cut association before optimization
            [~,vrtA] = computeCutAdjacency(cutIndices,objB,5);
            associatedCuts = findAssociatedCuts(triangleGroups,objCuts,current_centers,vrtA,cutCenters);
            
            % Optimize the cuts
            [objCuts,objCutNs,objCutKs,cutIndices] = optimizeCuts(cutIndices,objB,objN,objK,options);

            cutsIntersectBoundary = checkCutIntegrity(objCuts,objB);
            
            if cutsIntersectBoundary
                fprintf('cutIntersectsBoundary in object %d\n',obj)
            end
            
            
            
%                     figure(113)
%             clf(113)
% 
% %             for i = 1:size(objB,1)
% %                 if ~isnan(CreateObjectPartitionInfo.assignedCenter(i))
% %                     line([objB(i,1),current_centers(CreateObjectPartitionInfo.assignedCenter(i),1)], [objB(i,2),current_centers(CreateObjectPartitionInfo.assignedCenter(i),2)],'color',0.6*[1 1 1])
% %                 end
% %             end
%             
%             line(objB(:,1),objB(:,2),'Color','k','LineWidth',2)
%             hold on
%             plot(current_centers(:,1),current_centers(:,2),'ro','MarkerFaceColor','r')
%             plot(objCuts(:,1),objCuts(:,2),'o','LineWidth',1,'Color',0.3*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g')
%             plot(objCuts(:,3),objCuts(:,4),'d','LineWidth',1,'Color',0.3*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g')
% %             
% %             line(r_end{obj}(:,1)-topLeftB(1),r_end{obj}(:,2)-topLeftB(2),'Marker','*','Color','y','LineStyle','none')
% %             quiver(objB(:,1),objB(:,2),objN(:,1),objN(:,2),0.25,'Color',0.8*[1 1 1])
%             drawnow;
%             goDark(gcf)
%             for i = 1:size(objCuts,1)
%                 tmp = (objCuts(i,1:2) + objCuts(i,3:4))/2;
%                 line(objCuts(i,[1,3]),objCuts(i,[2,4]),'color',[0 0 1],'lineWidth',2)
%                 text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color','b')
%             end
%             
%             for i = 1:size(current_centers,1)
%                 text(current_centers(i,1),current_centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color','k')
%             end
%             
%             for i = 1:numel(CreateObjectPartitionInfo.partitions)
%                 line(CreateObjectPartitionInfo.partitions{i}(:,1),CreateObjectPartitionInfo.partitions{i}(:,2),'color','g')
%             end
% %             str = {'Note edge 4 should be connecting through the hole, but it is';
% %              'not because the normals at the hole and the boundary';
% %              'start/end points are perpendicular.'};
% %             title(str)
% %             grid on
% %             box on
% %             axis tight
%             daspect([1 1 1])
% %             error('somerror')
            
            
            
            
            
            % Get cut association after optimization
            [~,vrtA] = computeCutAdjacency(cutIndices,objB,5);
            associatedCutsAfterOptimization = findAssociatedCuts(triangleGroups,objCuts,current_centers,vrtA,cutCenters);
            
            % Combine associations before and after to get the complete
            % list.
            for i = 1:numel(associatedCuts)
                associatedCuts{i} = unique([associatedCuts{i};associatedCutsAfterOptimization{i}]);
            end
            
            % Create triangle cuts
            [objTriCuts,objTriCutNs,objTriCutKs,objTriCutInds,InfoTri] = createTriangleCuts(objB, objN, objK, triangleGroups, current_centers, associatedCuts, options);

            % Determine winning cuts
            [~,vrtA] = computeCutAdjacency(cutIndices,objB,1);            
            [objCuts,objCutNs,objCutKs,isTriCut] = chooseWinningCuts(objCuts,objCutNs,objCutKs,objTriCuts,objTriCutNs,objTriCutKs,associatedCuts,vrtA,objS,objI);
         
            % Remove any cuts that are not at concave points.
            if ~isempty(objCuts) && any(~isTriCut)
                toRemove = (mean(objCutKs,2) < 1/(2*options.maxRadius)) & ~isTriCut; % maybe use max instead of mean.
                objCuts(toRemove,:) = [];
            end
            
            % Create the set of pixels that the cuts cover. ----
            % First offset the cuts by the topLeftB
            objCuts(:,1:2) = objCuts(:,1:2) + topLeftB;
            objCuts(:,3:4) = objCuts(:,3:4) + topLeftB;


            % Get the linear indices of the pixels for each cut.
            current_cutInds = cell(size(objCuts,1),1);
            for current_cut = 1:size(objCuts,1)

                inds = rayTrace(objCuts(current_cut,1:2),objCuts(current_cut,3:4));
                current_cutInds{current_cut} = inds(:,1) + (inds(:,2)-1)*sizeI1;
            end

            cuts{obj} = cat(1,current_cutInds{:});
            Info{obj}.cutCalculationTime = toc(start1);
        else
            Info{obj}.cutCalculationTime = nan;
        end % if ~isempty(current_centers)

        catch ME
            toc(start)
            fprintf('Error in object %d\n',obj);
            rethrow(ME)
        end %try
    end %if need to declump
    
    % Feature extraction on each part of the object.
    
end % parfor
toc(start)
% profile off
% profile viewer
% reset(dev);
% gpuDevice([])

r0 = cat(1,r0{:});
r_end = cat(1,r_end{:});
BWtmp = BW;
BWtmp(cat(1,cuts{:})) = 0;

figure(2)
clf(2)
imshow(BWtmp)
hold on
plot(r0(:,2),r0(:,1),'r.')
plot(r_end(:,2),r_end(:,1),'b.')


% for i = 1:size(B{ind},1)
%     if ~isnan(assignedCenters(i))
%         line([B{ind}(i,2),centers{ind}(assignedCenters(i),2)], [B{ind}(i,1),centers{ind}(assignedCenters(i),1)],'color',0.6*[1 1 1])
%     end
% end
% line(B{ind}(isnan(assignedCenters),2),B{ind}(isnan(assignedCenters),1),'linestyle','none','marker','s','markerSize',5,'color','r')
centers = cat(1,centers{:});
plot(centers(:,2),centers(:,1),'bo','MarkerSize',3,'MarkerFaceColor','b')
% ylim([540 680])
% xlim([3420 3570])
drawnow;
goDark(gcf);


%%

% figure(4)
% clf(4)
% times = cellfun(@(x) x.solverTime,Info(~cellfun(@isempty,Info)));
% histogram(times)
% hold on
% times = cellfun(@(x) x.cutCalculationTime,Info(~cellfun(@isempty,Info)));
% histogram(times)
% xlabel('time')
% ylabel('counts')
% legend({'ode solver time','cut calculation time'})
% drawnow;
% goDark(gcf)
% for obj = 1:numel(B)
%     text(mean(B{obj}(:,2))+20,mean(B{obj}(:,1))+20,sprintf('%0.0f',10*abs(sum(curvature{obj}(curvature{obj}>0))/sum(abs(curvature{obj})))),'Color','g','Clipping','on','BackgroundColor',0.2*[1 1 1])
%     text(mean(B{obj}(:,2))+20,mean(B{obj}(:,1))+33,sprintf('%d',potentialParameters.extent(extentUsed(obj))),'Color','g','Clipping','on','BackgroundColor',0.2*[1 1 1])
% end




% p = hexplot([Info{1}.ComputeInitialPoints.CurvatureHexData.Y(:),Info{1}.ComputeInitialPoints.CurvatureHexData.X(:)],20,1,'maxHexSize',1,'Parent',gca,'orientation',90);
% p.FaceColor = 'none'; p.EdgeColor = 'r';

