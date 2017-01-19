function [cuts,dcuts,cutCntrs,cutIndices,partitions,Info] = calculateCuts(B,n,selectedCenter,centers)
% calculateCuts will compute the cuts to boundary B when each boundary
% element has been asigned to one of the N centers in selectedCenter.
%
% Input parameters:
% B - Mx2 array with the boundary coordinates
% n - Mx2 array with the boundary normals (pointing inward)
% selectedCenter - Mx1 array with the center that the boundary element has
%                  been asigned to
% centers - Nx2 array giving the centers
%
% Output parameters:
% cuts - Lx4 array where the i'th cut is between the vertices given by
%        cuts(i,1:2) and cuts(i,3:4)
% dcuts - Lx4 boundary normals at the cut vertices
% cutCntrs - Lx2 the center or centers associated with a cut
% cutIndices - Lx2 the indices of the cut vertices in the boundary B
% partions - 1xN cell array with the individual nuclei boundaries.
% Info - structure with fields
%        RebuildBoundary - information strucuture returned by
%                          rebiuldBoundary function
%        CleanCuts - information structure returned by cleanCuts function
%
%   Note, if there is only one input the, cuts, dcuts, cutCntrs,
%   cutIndices, and Info are all empty and partitions{1} = B.
%
% This function will break up the boundary into N partitions that are all
% clockwise orriented. The cuts necessary to break up the boundary B into
% these partitions are returned in cuts. Each cut is further labeled with
% the center or centers it separates.
%
% See also REBIULDBOUNDARY CLEANCUTS ORDERCUTVERTICES COMPUTECUTINDICES
% COMPUTECUTADJACENCY

% Note: A hard pair is two cuts such that vertex 1 of cut 1 is within
% HARD_PAIR_DISTANCE of vertex 1 of cut 2 and vertex 2 of cut 1 is within 
% HARD_PAIR_DISTANCE of vertex 2 of cut 2. The same goes for a soft pair,
% but the soft pair must not be a hard pair.

% James Kapaldo
% 2016-10-05


if nargout > 5
    debug = 1;
else
    debug = 0;
end

HARD_PAIR_DISTANCE = 1;
SOFT_PAIR_DISTANCE = 9;


N = size(centers,1);

% If there is only one center, then we have nothing to do.
if N == 1
    cuts = [];
    dcuts = [];
    cutCntrs = [];
    partitions{1} = B;
    Info = [];
    return;
end

% Initialized some arrays1
partitions = cell(1,N);
cuts = zeros(20,4);
dcuts = zeros(20,4);
cutCntrs = zeros(20,1);
counter = 1;

% Rebiuld each partition and get the cuts.
cols = {'r','b'};
for i = 1:N

    % Get the indices of the boundary assigned to the curernt center.
    currentAssignment = selectedCenter == i;
    partitionB = B(currentAssignment,:);
    partitionN = n(currentAssignment,:);
    
    % If the current set of boundary edges are cut by the start/end of the
    % boundary, then rotate the current boundary so that the only places we
    % have breaks are at the true breaks.
    if currentAssignment(1) && currentAssignment(end)

        % We know that the first and, since we are in this if
        % statement, last vertex of the boundary are in this partition.
        % We also know we still have a cut, which means
        % currentAssignment will have a form of
        %
        % currentAssignment = [1 1 ... 1 1 0 0 ... 0 0 1 1 ... 1 1]
        %
        % We are looking for the first 1 0 switch, so we can simply use
        % find(diff(~currentAssignment),1,'first')

        firstBreak = find(~currentAssignment,1,'first')-1;
        
        partitionB = circshift(partitionB,-firstBreak,1);
        partitionN = circshift(partitionN,-firstBreak,1);
    end
    
    if ~isempty(partitionB)
%         fprintf('partition %d\n',i)
        if debug
            [partitionB,partitionCuts,partitionCutsN,PartitionInfo] = rebuildBoundary(partitionB,partitionN,centers(i,:));
            % Initialize the debug structure if its the first iteration.
            if i == 1
                RebuildInfo(N) = PartitionInfo;  %#ok<AGROW>
            end
            RebuildInfo(i) = PartitionInfo; %#ok<AGROW>
        else
            [partitionB,partitionCuts,partitionCutsN] = rebuildBoundary(partitionB,partitionN,centers(i,:));
        end
        
        % It can be that the partition is completly close and there are not
        % cuts if the region is only connected with a 4-neighborhood to
        % another region. In this case we will find that the indices do
        % have a break, and we will make that our cut.
        if isempty(partitionCuts)
            partitionIdx = find(currentAssignment); % indices into B, not partitionB
            dPartitionIdx = circshift(partitionIdx,-1,1) - partitionIdx;
            cutLocationIdx = (dPartitionIdx ~= 1) & (dPartitionIdx ~= -(size(B,1)-1));
            cutLocationIdxp1 = find(circshift(cutLocationIdx,1,1));
            cutLocationIdx = find(cutLocationIdx);
            numberOfCuts = numel(cutLocationIdx);

            partitionCuts = zeros(numberOfCuts,4);
            partitionCutsN = zeros(numberOfCuts,4);

            for cutNum = 1:numel(cutLocationIdx)
                [partitionCuts(cutNum,:),partitionCutsN(cutNum,:)] = orderCutVertices(B(partitionIdx(cutLocationIdx(cutNum)),:),B(partitionIdx(cutLocationIdxp1(cutNum)),:),n(partitionIdx(cutLocationIdx(cutNum)),:),n(partitionIdx(cutLocationIdxp1(cutNum)),:));
            end
        end

        m = size(partitionCuts,1);
        
        if m+counter-1 > size(cuts,1)
            cuts = [cuts;zeros(20,4)]; %#ok<AGROW>
            dcuts = [dcuts;zeros(20,4)]; %#ok<AGROW>
            cutCntrs = [cutCntrs;zeros(20,1)]; %#ok<AGROW>
        end

        cuts(counter:m+counter-1,:) = partitionCuts;
        dcuts(counter:m+counter-1,:) = partitionCutsN;
        cutCntrs(counter:m+counter-1) = i;
        
        counter = counter + m;
        
        partitions{i} = [partitionB;partitionB(1,:)];
    end
end

cuts(counter:end,:) = [];
dcuts(counter:end,:) = [];
cutCntrs(counter:end) = []; % Note this gives the partition number.

% Create center association ==============================================
cutIndices = computeCutIndices(cuts,dcuts,B,n);
cutA = computeCutAdjacency(cutIndices,B,[HARD_PAIR_DISTANCE, SOFT_PAIR_DISTANCE]);
cutA = cellfun(@(x) triu(x)==2, cutA,'UniformOutput',0);

[hardPairs_i,hardPairs_j] = find(cutA{1});
[softPairs_i,softPairs_j] = find(cutA{2} & ~cutA{1});
hardPairs = [hardPairs_i,hardPairs_j];
softPairs = [softPairs_i,softPairs_j];

unPaired = find(~ismember((1:size(cuts,1))', [hardPairs(:); softPairs(:)]));

hardPairCntrs = cutCntrs(hardPairs);
softPairCntrs = cutCntrs(softPairs);
unPairedCntrs = [cutCntrs(unPaired,:), nan(numel(unPaired),1)];

if size(hardPairs,1)==1
    hardPairCntrs = hardPairCntrs';
end
if size(softPairs,1)==1
    softPairCntrs = softPairCntrs';
end

cutCntrAssociation = zeros(size(cuts,1),2);
for pair = 1:size(hardPairs,1)
    cutCntrAssociation(hardPairs(pair,:),:) = hardPairCntrs(pair,:) .* [1;1];
end
for pair = 1:size(softPairs,1)
    cutCntrAssociation(softPairs(pair,:),:) = softPairCntrs(pair,:) .* [1;1];
end
for pair = 1:size(unPaired,1)
    cutCntrAssociation(unPaired(pair),:) = unPairedCntrs(pair,:);
end

cutCntrs = cutCntrAssociation;

% figure(121)
% clf(121)
% line(B(:,1),B(:,2),'Color','k','LineWidth',2)
% hold on
% plot(centers(:,1),centers(:,2),'ro','MarkerFaceColor','r')
% 
% for i = 1:size(cuts,1)
%     tmp = (cuts(i,1:2) + cuts(i,3:4))/2;
%     line(cuts(i,[1,3]),cuts(i,[2,4]),'color',[0 0 1],'lineWidth',2)
%     text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color','b')
% end
% drawnow;
% goDark(gcf)


% Clean up the cuts ======================================================
options.HARD_PAIR_DISTANCE = HARD_PAIR_DISTANCE;
options.SOFT_PAIR_DISTANCE = SOFT_PAIR_DISTANCE;
[cuts,dcuts,cutCntrs,cutIndices,CleanCutsInfo] = cleanCuts(B,n,cuts,dcuts,cutCntrs,options);


% figure(122)
% clf(122)
% line(B(:,1),B(:,2),'Color','k','LineWidth',2)
% hold on
% plot(centers(:,1),centers(:,2),'ro','MarkerFaceColor','r')
% 
% for i = 1:size(cuts,1)
%     tmp = (cuts(i,1:2) + cuts(i,3:4))/2;
%     line(cuts(i,[1,3]),cuts(i,[2,4]),'color',[0 0 1],'lineWidth',2)
%     text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color','b')
% end
% drawnow;
% goDark(gcf)


% Save debug information if requested
if debug
    Info.CleanCuts = CleanCutsInfo;
    Info.RebuildBoundary = RebuildInfo;
end

end
% calculateCuts : changeLog
% 2016-10-14 : modified the code that gets the firstBreak
% 2016-10-18 : removed bug where a cut was used more than once as a hard
%              cut and a soft cut
% 2016-10-27 : major re-write of second half of function, created the
%              helper functions, doCutsIntersect, getIndicesToTheSide,
%              cleanCuts.


function answer = doCutsIntersect(cut1,cut2)
% DOCUTSINTERSECT  Determine if the two cuts intersect and return true or
% false.

% James Kapaldo
% 2016-10-27

    [~,~,t1,~] = intersections(cut1([1,3]),cut1([2,4]),cut2([1,3]),cut2([2,4]));    
    answer = ~isempty(t1);
end

function vrts = getIndicesToTheSide(idx,nanLocs)
% GETINDICESTOTHESIDE  Return the indices to either side of idx on the same
% contour that idx is on.
%
% This function assumes that the index comes from a boundary that can have
% several contours delimited by a row of nans.
%
% nanLocs should be the array [0, find(isnan(B(:,1))), size(B,1)+1]
% where B is the boundary array.

% James Kapaldo
% 2016-10-27

if numel(nanLocs)==2
    vrts = [idx-1,idx+1];
else
    endSide = nanLocs == idx+1;
    startSide = nanLocs == idx-1;

    if any(endSide) && any(startSide)
        % The contour has only one vertex. return empty
        vrts = [];
        return;
    end

    if any(endSide)
        % The index is the last vertex in a contour, need to wrap to
        % beginning.
        vrts = [idx-1,nanLocs(find(endSide)-1)+1];
        return;
    end

    if any(startSide)
        % The index is the first vertex in a contour, need to wrap to
        % end
        vrts = [nanLocs(find(startSide)+1)-1, idx+1];
        return;
    end
    vrts = [idx-1,idx+1];
end


end

function [cuts, cutNs, cutIndices] = tryToFixSelfAdjacency(B,n,cuts,cutNs,cutIndices,options)
% TRYTOFIXSELFADJACENCY
%
% It can be it two regions are only touching by a 4-connected
% neighborhood, that the outer boundary overlaps. In this case, the
% one (or two) cuts that are formed at this juncture can both be
% self-adjacent as one of the vertices could have two possible
% indices, a close one, and one on the opposite side of the
% contour. This function will check to see if the above situation
% occurs and fix it by use the index on the far side of the
% contour.

% James Kapaldo
% 2016-10-30

% Find any contour delimiters
nanLoc = find(isnan(B(:,1)));
if isempty(nanLoc)
    nanLoc = size(B,1)+1;
end

% Get the unique vertices of the outer contour.
[Bu,~,ic] = unique(B(1:nanLoc(1)-1,:),'rows');

% The case above is only true if there are repeated vertices in the
% outer contour.
if (nanLoc(1)-1) == size(Bu,1)
    return;
end

cutA = computeCutAdjacency(cutIndices,B,options.SOFT_PAIR_DISTANCE);
selfAdjacent = logical(diag(cutA));

saCuts = cuts(selfAdjacent,:);
saCutInds = cutIndices(selfAdjacent,:);

% Get the duplicate vertices
vrts = Bu(accumarray(ic,1)>1,:);

% Get the indices of these vertices in the original boundary
[vrtsIdx,~] = find( permute(all((B - permute(vrts,[3,2,1]))==0,2), [1 3 2]) );

vrts = B(vrtsIdx,:);

% Replace each duplicate vertex index with the other possible index
% (we can always replace since we already know the cut is
% self-adjacent)
for i = 1:size(saCuts,1)
    isvrt1 = sum((saCuts(i,1:2) - vrts)==0,2)==2;
    isvrt2 = sum((saCuts(i,3:4) - vrts)==0,2)==2;
    
    if any(isvrt1) && any(isvrt2)
        % we only switch one of them
        switched = 0;
        
        vrt1idxs = vrtsIdx(isvrt1,:); % this should be at most 2x1
        newInd = vrt1idxs(vrt1idxs~=saCutInds(i,1));
        if ~any(cutIndices(:)==newInd)
            saCutInds(i,1) = newInd;
            switched = 1;
        end
        
        if ~switched
            vrt2idxs = vrtsIdx(isvrt2,:); % this should be at most 2x1
            newInd = vrt2idxs(vrt2idxs~=saCutInds(i,2));
            if ~any(cutIndices(:)==newInd)
                saCutInds(i,2) = newInd;
            end
        end
    else
    
        if any(isvrt1)
            vrt1idxs = vrtsIdx(isvrt1,:); % this should be at most 2x1
            newInd = vrt1idxs(vrt1idxs~=saCutInds(i,1));            
            if ~any(cutIndices(:)==newInd)
                saCutInds(i,1) = newInd;
            end
        end

        if any(isvrt2)
            vrt2idxs = vrtsIdx(isvrt2,:); % this should be at most 2x1
            newInd = vrt2idxs(vrt2idxs~=saCutInds(i,2));
            if ~any(cutIndices(:)==newInd)
                saCutInds(i,2) = newInd;
            end
        end
    end
end


saCuts = [B(saCutInds(:,1),:),B(saCutInds(:,2),:)];
saCutNs = [n(saCutInds(:,1),:),n(saCutInds(:,2),:)];

cuts(selfAdjacent,:) = saCuts;
cutNs(selfAdjacent,:) = saCutNs;
cutIndices(selfAdjacent,:) = saCutInds;

end


function [cuts,cutNs,cutCenters,cutIndices,Info] = cleanCuts(B,n,cuts,cutNs,cutCenters,options)
% CLEANCUTS  Ensure no vertices are shared and remove cuts whos vertices
% are adjacent to each other. Also, remove a cut from a pair of hard cuts.
%
% If two cuts share a vertex and one of the cuts is paired with another
% cut, then remove the paired cut. If both cuts are unpaired, then try to
% move the duplicate vertex over one, if we cannot move it over because
% that would create another duplicate vertex, then remove one of the cuts.
%
% Note: A hard pair is two cuts such that vertex 1 of cut 1 is within
% HARD_PAIR_DISTANCE of vertex 1 of cut 2 and vertex 2 of cut 1 is within 
% HARD_PAIR_DISTANCE of vertex 2 of cut 2. The same goes for a soft pair,
% but the soft pair must not be a hard pair.

% James Kapaldo
% 2016-10-27
Info.unpairedCutRemoved = false;

% Find the contour breaks
nanLocs = [0;find(isnan(B(:,1)));size(B,1)+1];

% Remove a duplicate cuts
[cuts,ia] = unique(cuts,'rows');
cutNs = cutNs(ia,:);
cutCenters = cutCenters(ia,:);

% Get the indices and adjacency
cutIndices = computeCutIndices(cuts,cutNs,B,n);

if size(cuts,1)==1
    return;
end

[cuts, cutNs, cutIndices] = tryToFixSelfAdjacency(B,n,cuts,cutNs,cutIndices,options);

[cutA,vrtA] = computeCutAdjacency(cutIndices,B,[0, options.HARD_PAIR_DISTANCE, options.SOFT_PAIR_DISTANCE]);


% Short self adjacent cuts
selfAdjacent = find(diag(cutA{3}));

% Get shared, hard, and soft info
cutA_share = triu(cutA{1});
vrtA_share = triu(vrtA{1});

cutA_hard = cutA{2}==2;
cutA_soft = cutA{3}==2 & ~cutA_hard;

cutA_share(selfAdjacent,:) = [];
vrtA_share(selfAdjacent,:) = [];
cutA_hard(selfAdjacent,:) = [];
cutA_soft(selfAdjacent,:) = [];

cutA_share(:,selfAdjacent) = [];
vrtA_share(:,selfAdjacent) = [];
cutA_hard(:,selfAdjacent) = [];
cutA_soft(:,selfAdjacent) = [];

cuts(selfAdjacent,:) = [];
cutNs(selfAdjacent,:) = [];
cutCenters(selfAdjacent,:) = [];
cutIndices(selfAdjacent,:) = [];

% Reshape indices to row (to coorespond with vertex indices)
cutIndicesRow = reshape(cutIndices',1,2*size(cuts,1));

[share_i,share_j] = find(cutA_share);
toRemove = false(size(cuts,1),1);

for vrtNum = 1:numel(share_i)
    
    associatedCutTypes = [any(cutA_hard(share_i(vrtNum),:)), any(cutA_hard(share_j(vrtNum),:)); 
                          any(cutA_soft(share_i(vrtNum),:)), any(cutA_soft(share_j(vrtNum),:))];
    
    if any(associatedCutTypes(1,:))
        % at least one cut has a hard pair, remove it
        if associatedCutTypes(1,1)
            toRemove(share_i(vrtNum)) = true;
        else
            toRemove(share_j(vrtNum)) = true;
        end
    elseif any(associatedCutTypes(2,:))
        % at least one cut has a soft pair. if both are soft pairs, remove
        % the cut with the closer soft pair. if one is soft and one is
        % unpaired, remove the soft pair
        
        if all(associatedCutTypes(2,:)) % both soft           
            Di = sum((cuts(share_i(vrtNum),:) - cuts(cutA_soft(share_i(vrtNum),:),:)).^2,2);
            Dj = sum((cuts(share_j(vrtNum),:) - cuts(cutA_soft(share_j(vrtNum),:),:)).^2,2);
            if Dj < Di
                toRemove(share_i(vrtNum)) = true;
            else
                toRemove(share_j(vrtNum)) = true;
            end
        else
            if associatedCutTypes(2,1)
                toRemove(share_i(vrtNum)) = true;
            else
                toRemove(share_j(vrtNum)) = true;
            end
        end
    else
        % both cuts are unpaired.
        vidx_i = 2*share_i(vrtNum)+[-1,0];
        vidx_j = 2*share_j(vrtNum)+[-1,0];
        
        sharedVrtIdx = cutIndicesRow(vidx_i(any(vrtA_share(vidx_i,vidx_j),2))); % This line assumes that two cuts do not lie directly on top of each other -- there can only be one shared vertex
        
        sideVrtIdx = getIndicesToTheSide(sharedVrtIdx,nanLocs);
        
        sideVrtsAreCuts = ismember(sideVrtIdx,cutIndicesRow);
        
        if all(sideVrtsAreCuts)
            % We cannot move a cut vertex over, so we will just throw away
            % one of the unpaired cuts.
            toRemove(share_i(vrtNum)) = true;
            Info.unpairedCutRemoved = true;
        else
            % One of the sides is open
            sideVrtIdx = sideVrtIdx(find(~sideVrtsAreCuts,1,'first'));
            cut1 = cuts(share_i(vrtNum),:);
            cut2 = cuts(share_j(vrtNum),:);
            
            tmpCut = cut1;
            tmpCut((find(cutIndices(share_i(vrtNum),:)==sharedVrtIdx)-1)*2 + (1:2)) = B(sideVrtIdx,:);
            
            if doCutsIntersect(tmpCut,cut2)
                tmpCut = cut2;
                tmpCut((find(cutIndices(share_j(vrtNum))==sharedVrtIdx)-1)*2 + (1:2)) = B(sideVrtIdx,:);
                if doCutsIntersect(tmpCut,cut1)
                    error('cuts intersect no matter which point is used, this should not happen, there is likely a programming but somewhere')
                end
                cuts(share_j(vrtNum),(find(cutIndices(share_j(vrtNum),:)==sharedVrtIdx)-1)*2 + (1:2)) = B(sideVrtIdx,:);
                cutNs(share_j(vrtNum),(find(cutIndices(share_j(vrtNum),:)==sharedVrtIdx)-1)*2 + (1:2)) = n(sideVrtIdx,:);
                cutIndices(share_j(vrtNum),cutIndices(share_j(vrtNum),:)==sharedVrtIdx) = sideVrtIdx;
            else
                cuts(share_i(vrtNum),(find(cutIndices(share_i(vrtNum),:)==sharedVrtIdx)-1)*2 + (1:2)) = B(sideVrtIdx,:);
                cutNs(share_i(vrtNum),(find(cutIndices(share_i(vrtNum),:)==sharedVrtIdx)-1)*2 + (1:2)) = n(sideVrtIdx,:);
                cutIndices(share_i(vrtNum),cutIndices(share_i(vrtNum),:)==sharedVrtIdx) = sideVrtIdx;
            end
        end
    end
end

% Remove hard pairs
cutA = computeCutAdjacency(cutIndices,B,options.HARD_PAIR_DISTANCE);
cutA = triu(cutA)==2;

toRemove = toRemove | any(cutA,2);

cuts(toRemove,:) = [];
cutNs(toRemove,:) = [];
cutIndices(toRemove,:) = [];
cutCenters(toRemove,:) = [];

end
% cleanCuts : changeLog
% 2016-10-30: added in fix for self-adjacent cuts where the vertices are
%             duplicates in the boundary (added function
%             tryToFixSelfAdjacency)



















% % % 
% % % 
% % % 
% % % % figure(113)
% % % % clf(113)
% % % % for i = 1:N
% % % %     line(partitions{i}(:,1),partitions{i}(:,2),'color','g')
% % % % end
% % % % hold on
% % % % plot(centers(:,1),centers(:,2),'ro')
% % % % plot(cuts(:,1),cuts(:,2),'ys')
% % % % plot(cuts(:,3),cuts(:,4),'c*')
% % % % drawnow;
% % % % goDark(gcf)
% % % 
% % % 
% % % % It could be that two cuts share the same vertex. This will lead to errors
% % % % when each cut is optimized - each vertex must be unique. 
% % % %
% % % % Determine if there are any duplicate vertices and determine what cuts
% % % % they are in.
% % % 
% % % % reshape the cuts as Nx4 -> 2*Nx2
% % % cuts = reshape(cuts',2,2*(counter-1))'; 
% % % 
% % % % Find the shared vertices.
% % % [~,~,ic1] = unique(cuts,'rows');
% % % doubleVrts = accumarray(ic1,1)>1;
% % % sharedVrts = sum((ic1 == find(doubleVrts')) .* (1:sum(doubleVrts)), 2); % array of size 3*Q x 1, each pair will have a different number.
% % % 
% % % % reshape them back to original form
% % % cuts = reshape(cuts',4,counter-1)'; % Nx4
% % % sharedVrts = reshape(sharedVrts',2,counter-1)'; % Nx2
% % % 
% % % 
% % % % Find the cuts that that are paired with each other. Two cuts are hard
% % % % pairs with each other if both of their end points are within
% % % % sqrt(HARD_PAIR_DISTANCE) of each other. Two cuts are soft pairs with each
% % % % other if there end points are within sqrt(SOFT_PAIR_DISTANCE). Hard pair
% % % % cuts will be replaced with a single cut with the shortest distance.
% % % %
% % % % These pairs are not just to replace hard pairs with one cut, but are also
% % % % for labeling the centers that a cut seperates. Cut pairs (hard or soft)
% % % % are associated with two centers (the center on either side of the cut).
% % % % If a cut is not in a hard or soft pair, then it will be called unpaired,
% % % % and only one center is associated with it. (Keeping track of this will
% % % % be important when triangle optimization is done later on and we need to
% % % % find all edges associated with a set of centers.)
% % % 
% % % pdistInds = getPdistInds(counter-1);
% % % 
% % % pairDists = pdist(cuts)';
% % % 
% % % % Get indices of hard pairs
% % % hardPairInds = pairDists <= HARD_PAIR_DISTANCE;
% % % 
% % % % Get indices of soft pairs and make sure that no cut that has already
% % % % been declared as a hard pair is repeated.
% % % softPairInds = (pairDists >  HARD_PAIR_DISTANCE) & (pairDists <= SOFT_PAIR_DISTANCE) & (~any(any(pdistInds == reshape(pdistInds(hardPairInds,:),1,1,[]),3),2));
% % % 
% % % % Hard-, soft-, and un-paired indices
% % % hardPairs = pdistInds(hardPairInds,:);
% % % softPairs = pdistInds(softPairInds,:);
% % % unPaired = find(~ismember((1:counter-1)', [hardPairs(:); softPairs(:)]));
% % % 
% % % 
% % % hardPairCntrs = cutCntrs(hardPairs);
% % % softPairCntrs = cutCntrs(softPairs);
% % % unPairedCntrs = [cutCntrs(unPaired,:), nan(numel(unPaired),1)];
% % % 
% % % if size(hardPairs,1)==1
% % %     hardPairCntrs = hardPairCntrs';
% % % end
% % % if size(softPairs,1)==1
% % %     softPairCntrs = softPairCntrs';
% % % end
% % % 
% % % 
% % % cutCntrAssociation = zeros(size(cuts,1),2);
% % % for pair = 1:size(hardPairs,1)
% % %     cutCntrAssociation(hardPairs(pair,:),:) = hardPairCntrs(pair,:) .* [1;1];
% % % end
% % % for pair = 1:size(softPairs,1)
% % %     cutCntrAssociation(softPairs(pair,:),:) = softPairCntrs(pair,:) .* [1;1];
% % % end
% % % for pair = 1:size(unPaired,1)
% % %     cutCntrAssociation(unPaired(pair),:) = unPairedCntrs(pair,:);
% % % end
% % % 
% % % 
% % % pairData = [hardPairs, ones(size(hardPairs,1),1); 
% % %             softPairs, 2*ones(size(softPairs,1),1);
% % %             unPaired, unPaired, 3*ones(size(unPaired,1),1)];
% % % 
% % % toRemove = false(size(cuts,1),1);
% % % for vrt = 1:max(sharedVrts(:))
% % %     cutNums = find(any(sharedVrts==vrt,2)); % 2x1
% % %     pairIdx = find(any(any(pairData(:,1:2) == permute(cutNums,[3,2,1]),3),2)); % size(pairData,1)x1
% % %     pairTypes = pairData(pairIdx,3) == permute(1:3,[1,3,2]); % size(pairData,1)x1x3
% % %     
% % %     if any(pairTypes(:,1,1)) % there is at least one hard pair
% % %         % Either both are hard pairs, or one is a hard pair and the other
% % %         % is either a soft pair or unpaired.
% % %         % Need to remove one of the hard cuts.
% % %         
% % %         if all(pairTypes(:,1,1))
% % %             toRemove(cutNums(1)) = true;
% % %         else
% % %             % We want the pairData of the pairIdx that is a hard pair. this
% % %             % will give the 2 cuts that make up the pair. compair these two
% % %             % cuts with the cutNums (which are the cuts with the duplicate
% % %             % vertex) to find the cut with the duplicate vertex. set it to
% % %             % removal
% % % 
% % %             toRemove(cutNums(any(cutNums == pairData(pairIdx(pairTypes(:,1,1)),1:2),2))) = true; % start at the very beginning..(or inside as it were)
% % %         end
% % %         
% % %     elseif any(pairTypes(:,1,2)) % there is at least one soft pair
% % %         % Either both cuts are soft pairs or one is a soft pair and one is
% % %         % unpaired.
% % %         % Need to remove one of the soft pairs.
% % %         
% % %         if all(pairTypes(:,1,2))
% % %             % remove the longer cut
% % %             d1 = sum((cuts(cutNums(1),1:2) - cuts(cutNums(1),3:4)).^2,2);
% % %             d2 = sum((cuts(cutNums(2),1:2) - cuts(cutNums(2),3:4)).^2,2);
% % %             toRemove(cutNums((d1<d2)+1)) = true;
% % %         else
% % %             % We want the pairData of the pairIdx that is a soft pair. this
% % %             % will give the 2 cuts that make up the pair. compare these two
% % %             % cuts with the cutNums (which are the cuts with the duplicate
% % %             % vertex) to find the cut with the duplicate vertex. set it to
% % %             % removal
% % %             toRemove(cutNums(any(cutNums == pairData(pairIdx(pairTypes(:,1,2)),1:2),2))) = true; % start at the very beginning..(or inside as it were)
% % %         end
% % %         
% % %     else % both cuts are unapred
% % %         % This case should be very rare, and properly fixing this case
% % %         % would be quite .. in depth .. An outline for how one might fix it
% % %         % is given below, but for now (2016-10-18) it will not be
% % %         % implemented. One of the cuts will simply be removed and a warning
% % %         % message displayed.
% % %         
% % %         % Find all indices of the cut vertices - computeCutIndices()
% % %         % Determine if there is at least one index to the side of the
% % %         %   repeaded index that does not have a cut vertex. If there is a
% % %         %   free index, then find the unpaired cut on that side and move
% % %         %   the duplicate vertex to the free vertex. (Remember to keep in
% % %         %   mind that there could be multiple contours seperated by NaN's
% % %         %   so proper circular shifting of each individual contour must be
% % %         %   taken into account.)
% % %         % If both indices to either side of the duplicate index are already
% % %         %   full, then you will need to move one of the occupied indices
% % %         %   over to create a free index next to the duplicate index. With
% % %         %   a free index now available, follow the step just above.
% % %         
% % %         toRemove(cutNums(1)) = true;
% % % %         warning('calculateCuts:duplicateVertexInUnpairedCuts','Two unpaired cuts share the same vertex. One of the cuts will be removed, if this is not what you want, then you will have to implement proper handeling of this error yourself.')
% % %         
% % %     end
% % % end
% % % cuts(toRemove,:) = [];
% % % dcuts(toRemove,:) = [];
% % % cutCntrAssociation(toRemove,:) = [];
% % % 
% % % 
% % % % Get the hard pairs again (this is simpler than tracking them above, but
% % % % will take a very small amount of time more)
% % % %
% % % % Now we could replace each hard pair with a single cut with the minimum
% % % % distance between the 4 vertices; however, all cuts will be optimized
% % % % after this function and the optimziation search radius should be at least
% % % % 1 pixel; therefore, the cut will move to the best position anyway.
% % % % Therefore, we will simply delete one cut of each hard pair.
% % % 
% % % if size(cuts,1)>1
% % %     pdistInds = getPdistInds(size(cuts,1));
% % % 
% % %     pairDists = pdist(cuts)';
% % % 
% % %     % Get indices of hard pairs
% % %     hardPairs = pdistInds(pairDists <= HARD_PAIR_DISTANCE,1);
% % %     cuts(hardPairs,:) = [];
% % %     dcuts(hardPairs,:) = [];
% % %     cutCntrAssociation(hardPairs,:) = [];
% % % end
% % % 
% % % cutCntrs = cutCntrAssociation;

% [hardPairInds, (pairDists >  HARD_PAIR_DISTANCE) & (pairDists <= SOFT_PAIR_DISTANCE), softPairInds, pdistInds]
% 
% 
% 
% N_hardPairs = size(hardPairs,1);
% 
% 
% if N_hardPairs==1
%     hardPairCntrs = hardPairCntrs';
% end
% 
% newCuts = zeros(N_hardPairs,4);
% newdCuts = zeros(N_hardPairs,4);
% 
% for i = 1:N_hardPairs
%     points = cuts(hardPairs(i,:),:);
%     dpoints = dcuts(hardPairs(i,:),:);
%     % points = [v1 v2;
%     %           v3 v4];
%     % point v1 and v3 are within ~sqrt(HARD_PAIR_DISTANCE) of each other
%     % and v2 and v4 are withing ~sqrt(HARD_PAIR_DISTANCE) of each other.
%     dists = sum((points(:,1:2) - permute(points(:,3:4),[3,2,1])).^2,2);
%     % dists = (page 1) [d12; d32] (page 2) [d14 d34]
%     [~,minDist] = min(dists(:));
%     switch minDist
%         case 1
%             newCuts(i,:) = points(1,:);
%             newdCuts(i,:) = dpoints(1,:);
%         case 2
%             newCuts(i,:) = [points(2,1:2), points(1,3:4)];
%             newdCuts(i,:) = [dpoints(2,1:2), dpoints(1,3:4)];
%         case 3
%             newCuts(i,:) = [points(1,1:2), points(2,3:4)];
%             newdCuts(i,:) = [dpoints(1,1:2), dpoints(2,3:4)];
%         case 4
%             newCuts(i,:) = points(2,:);
%             newdCuts(i,:) = dpoints(2,:);
%     end
% end