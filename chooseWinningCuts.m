function [cuts,cutNs,cutKs,isTriangleCut,Scores] = chooseWinningCuts(cuts,cutNs,cutKs,triCuts,triCutNs,triCutKs,associatedCuts,vrtA,S,I)
% CHOOSEWINNINGCUTS  Decide between the triangle cuts and the cuts
% associated with the triangles
%
% [cuts,cutNs,cutKs,isTriangleCut,Scores] = chooseWinningCuts(cuts,cutNs,cutKs,triCuts,triCutNs,triCutKs,associatedCuts,vrtA,S,I)
% 
% Input paramters : 
% cuts - Lx4 array where the i'th cut is between the vertices cuts(i,1:2)
%        and cuts(i,3:4)
%
% cutNs - Lx4 array giving the boundary normals at the cut vertices
%
% cutKs - Lx2 array giving the boundary curvature at the cut vertices
%
% triCuts - cell array with cuts for each triangle group
%
% triCutNs - cell array with boundary normals for each cut vertex in
%            triCuts
%
% triCutKs - cell array with boundary curvatures for at cut vertex in
%            triCuts
%
% associatedCuts - cell array with the indices of the cut vertices
%                  associated with each triangle group
%
% vrtA - vertex association matrix
%
% S - image edges
%
% I - 1 over the image intensity
%
% Output parameters :
% cuts - final set of cuts
% cutNs - final set of normals for cuts
% cutKs - final set of curvatures for cuts
% isTriangleCut - logical array. The i'th element is true if cuts(i,:) came
%                 from a triangle cut
% Scores - cell array of the scores used in determining if the triangle
%          cuts should be used or not (first column is triangle cut scores)
%
% The dot product between the cut directions and boundary normals at the
% cut vertices, the boundary curvature at the cut vertices, the mean edge
% intensity along the cuts, and the mean image intensity (acutally 1/image
% intensity) along the cuts are calculated for each group. If the triangle
% cuts have three of these values larger than the normal cuts, then the
% triangle cuts are used. If its a tie and the triangle cuts with two of
% these and the normal cuts win two of these, then whichever one has the
% higher total score will be chosen.
%
% See also CREATEOBJECTPARTITIONS CREATETRIANGLECUTS

% James Kapaldo
% 2016-10-19

if isempty(triCuts)
    isTriangleCut = false(size(cuts,1),1);
    Scores = [];
    return;
end


imSize1 = size(S,1);

winningCuts = cell(1,numel(triCuts));
winningCutNs = cell(1,numel(triCuts));
winningCutKs = cell(1,numel(triCuts));

cutsToRemove = false(size(cuts,1),1);

vrts = reshape(1:2*size(cuts,1), 2, size(cuts,1))';

Scores = cell(1,numel(triCuts));

for i = 1:numel(triCuts)

    triCutsToConsider = ~isnan(triCutNs{i}(:,1));
    cutsToConsider = associatedCuts{i};    
    vrtsToConsider = vrts(cutsToConsider,:);
    
    % Triangle cut scores -----------------------------------------------
    rt = triCuts{i}(triCutsToConsider,3:4) - triCuts{i}(triCutsToConsider,1:2); %Lx2
    nt = triCutNs{i}(triCutsToConsider,1:2);
    kt = triCutKs{i}(triCutsToConsider,1);
    
    dt = sqrt(sum(rt.^2,2));
    rt = rt ./ dt;
    at = sum(rt .* nt,2);

%     triD = mean(dt)*sqrt(3);
    triA = mean(at);
    triK = mean(kt);
    
    % Image edge overlap
    triCutPix = cell(size(triCuts{i},1),1);
    for cutNum = 1:size(triCuts{i},1)
        inds = rayTrace(triCuts{i}(cutNum,1:2),triCuts{i}(cutNum,3:4));
        triCutPix{cutNum} = inds(:,1) + (inds(:,2)-1)*imSize1;
    end

    triO = mean(S(cat(1,triCutPix{:})));
    triIO = mean(I(cat(1,triCutPix{:})));
    
    % Cut scores --------------------------------------------------------
    r1 = cuts(cutsToConsider,1:2);
    n1 = cutNs(cutsToConsider,1:2);
    
    r2 = cuts(cutsToConsider,3:4);
    n2 = cutNs(cutsToConsider,3:4);
    
%     k = objCutKs(cutsToConsider,:);
    
    r = r1 - r2;
    d = sqrt(sum(r.^2,2));
    r = r ./ d;
    
    a1 = sum(-r .* n1,2);
    a2 = sum( r .* n2,2);
    
%     D = mean(d);
    A = mean([a1;a2]);
    K = computeMeanCurvature(reshape(cutKs',numel(cutKs),1),vrtsToConsider,vrtA);
%     K = mean(k(:));
    
    % Image edge overlap
    cutPix = cell(numel(cutsToConsider),1);
    for cutNum = 1:numel(cutsToConsider)
        inds = rayTrace(cuts(cutsToConsider(cutNum),1:2),cuts(cutsToConsider(cutNum),3:4));
        cutPix{cutNum} = inds(:,1) + (inds(:,2)-1)*imSize1;
    end
    
    O = mean(S(cat(1,cutPix{:})));
    IO = mean(I(cat(1,cutPix{:})));
    % Normalized old and new scores (by type) by their mean
    scores = [triO, O;
              triIO, IO;
              triA, A; 
              triK, K];
    scores = round( scores ./ mean(scores,2), 3);
    Scores{i} = scores;
%     scores = round(scores,3)
%     score = sum(scores,1)
    numTriWinners = sum(scores(:,1) > scores(:,2));
    
    if numTriWinners >= 3
        cutsToRemove(cutsToConsider) = true;
        winningCuts{i} = triCuts{i};
        winningCutNs{i} = triCutNs{i};
        winningCutKs{i} = triCutKs{i};
    elseif numTriWinners == 2
        score = sum(scores,1);
        if score(1) > score(2)
            cutsToRemove(cutsToConsider) = true;
            winningCuts{i} = triCuts{i};
            winningCutNs{i} = triCutNs{i};
            winningCutKs{i} = triCutKs{i};
        end
    end
    
end

 
cuts = [cat(1,winningCuts{:}); cuts(~cutsToRemove,:)];
cutNs = [cat(1,winningCutNs{:}); cutNs(~cutsToRemove,:)];
cutKs = [cat(1,winningCutKs{:}); cutKs(~cutsToRemove,:)];

isTriangleCut = true(size(cuts,1),1);
isTriangleCut(1:sum(~cutsToRemove)) = false;
isTriangleCut = flip(isTriangleCut);



end


function K = computeMeanCurvature(k,vrtInds,vrtA)
% COMPUTEMEANCURVATURE  Compute the mean curvature taking into account
% adjacent vertices.
%
% The two curvature values of adjacent vertices will be replace with a
% single value, the mean of the adjacent vertices' curvature. 
%
% Why do this? Consider the normal cuts, the likely form a triangle such
% that the 6 cut vertices are actually only in 3 locations. The curvature
% in these 3 locations is being double counted. Then note that the triangle
% cuts will only have one vertex at these three locations. Thus, replacing
% the curvature of adjacent vertices with their mean lets us better compair
% the curvature between the normal and triangle cuts.

ks = zeros(numel(vrtInds),1);
vrtInds = vrtInds(:);

% fprintf('============================\n')
% fprintf('computeMeanCurvature:\n')
% k(vrtInds)
% mean(k(vrtInds))
counter = 1;
while ~isempty(vrtInds)
%     fprintf('counter = %d\n',counter)
    v = vrtInds(1);
    allAdjVs = find(vrtA(v,:));
    adjVs = any(vrtInds == allAdjVs,2);
%     v
%     vrtInds(adjVs)
%     k([v;vrtInds(adjVs)])
    ks(counter) = mean(k([v;vrtInds(adjVs)]));
    
    vrtInds([1;find(adjVs)]) = [];
    counter = counter + 1;
end

ks(counter:end) = [];
% ks
K = mean(ks);

% fprintf('End : computeMeanCurvature\n')
% fprintf('============================\n')


end