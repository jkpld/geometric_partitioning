function [B,cuts,dcuts,Info] = rebuildBoundary(B,n,center)
% REBUILDBOUNDARY  Take in the pieces of boundary assigned to a center and
% create a well oriented contour.
%
% Assume you had good boundary. Now assume you took N of the elements from
% somewhere in the center and moved them to a different location. Example:
% let the below be an index array
% 
% good array       array with shifted group
%  1                1
%  2                2
%  3                3
%  4                8
%  5                9
%  6               10
%  7                4
%  8                5
%  9                6
% 10                7
%
% This function will take in the array with the shifted group and rebuild
% it to be the good array.
% This function will assume that each peice of the input array has the same
% orrientation. (That is that the shifted group in the example above is NOT
% [7;6;5;4].)
%
% B is the input boundary and should have two columns (x | y).
% dB is the normals (pointing inwards) of the boundary (dx | dy).
%
% The input boundary should not be "closed" --- the first and last
% point should not be the same.
%
% This function also assumes that the boundary is made on an image -- that
% is that the distance between connected points should have a distance <=
% sqrt(2).
%
% cuts will be an array of size n x 4, where n is the number of cuts (line
% segments longer than sqrt(2)). The columns will have the form
% [x1 y1 x2 y2] where [x1 y1] is one end of the cut and [x2 y2] is the
% other end of the cut. If the cut is more horizontal, then x1 <= x2; if
% the cut is more vertical, then y1 <= y2.
%
% dcuts will be the boundary normals at the cut locations.
%
% Info - a structure array with the fields
%       originalOrder - original order of the boundary pieces
%       bounds - boundaries of the boundary pieces
%       starts - start vertices of the boundary pieces
%       ends - end vertices of the boundary pieces
%       newOrder - order of the boundary pieces after ordering them
%       finalOrder - order of the boundary pieces after taking into accound
%                    far-away peices and double pieces.
%
% See also ORDERCUTVERTICES

% James Kapaldo
% 2016-10-05

SAME_DIRECTION_THRESHOLD = 0.9;

if nargout > 3
    debug = 1;
    Info = struct('originalOrder',[],'bounds',[],'starts',[],'ends',[],'finalOrder',[],'newOrder',[]);
else
    debug = 0;
end

cuts = [];
dcuts = [];

% Get the square distance between adjacent boundary vertices..
d = sum(diff(B,[],1).^2,2);

% Find the groups of edges that make up the polygon by looking for any two
% adjacet boundary vertices with a distance larger than 2.
groups = find(d > 2);

% If there are not groups, then our boundary only had one group. We just
% need to see if the start point and end point are also right next to each
% other, and, if not, then add in one cut from the start to the end.
if isempty(groups)
    if sum((B(end,:) - B(1,:)).^2) > 2
        [cuts,dcuts] = orderCutVertices(B(1,:),B(end,:),n(1,:),n(end,:));
    end
    return;
end

% Get the boundaries of the groups.
bnds = [0; groups; size(B,1)];

% Get the start and end points of each group.
starts = B([1;groups+1],:);
ends = B([groups;end],:);

ds = n([1;groups+1],:);
de = n([groups;end],:);

% Iterate through the number of groups. Start with the first group and find
% the next group whose start point is closest to the current groups end
% point. Do this the-number-of-groups + 1 times. It could happen that one
% group of the boundary is very far from all other groups. If this group
% that is far away is not the first group, then it would never be selected
% as the next group. If this group was the first element then we would
% never come back to this element. We do not want to include these far away
% groups because they are likely just from a bad boundary-to-center
% assignement caused by a maximum radius distance that is too large. By
% iterating through the-number-of-groups + 1 times, if the far away group
% was first, then we will know it is far away since we do not come back to
% it on the last iteration. If the far away group was not first, then we
% will never get to it and we will repeat groups.

% -----------------------------------
% Don't just find the closest group, but find the group who maximizes the
% dot product between the boundary normal and the cut and minimizes the
% distance. This will be done by dividing the dot product by the distance,
% and then maximizing.
% -----------------------------------

% Start with the group whose start location is closest to the center.
[~,currentGroup] = min(sum((center-starts).^2,2));

% Initialize the current order of the groups and create an array to store
% the new order.
originalOrder = 1:numel(groups)+1;
originalOrder = circshift(originalOrder,-currentGroup+1);
newOrder = originalOrder;

% Save some information about the rebuilding for debugging purposes.
if debug
    Info.originalOrder = originalOrder;
    Info.bounds = bnds;
    Info.starts = starts;
    Info.ends = ends;
end

for i = 1:numel(groups)+1
%     fprintf('======================\n')
%     fprintf('group %d\n', i)
    % Determine if we have already been to the current group. If we have
    % been to the current group, then we are currently in a loop and there
    % is no reason to continue.
    if any(currentGroup == newOrder(1:i-1))
        % We do not remove the duplicate group point. This way, the first
        % and last group should always be the same. 
        newOrder(i+1:end) = []; 
        break;
    end
    
    % Get the distance from the end of the current group to the starts of
    % all other groups.
    
    r = ends(currentGroup,:) - starts; % Lx2 : L = size(starts,1)
    d = sqrt(sum(r.^2,2)); % Lx1
    
    r = r./d; % Lx2
    
    % Get the dot products between each normal vector and the cut
    % vector.

    a1 = sum(r.*ds,2); % Lx1
    a2 = sum(-r .* de(currentGroup,:),2); % Lx1
    
    % Maximize the value of the dot product divided by the distance.
    
    valueToOptimize = (a1+a2)./d; % Lx1 : I want to maximise the numerator and minimize the denominator
    [sortedValues,ordD] = sort(valueToOptimize,'descend');
    
    % It can happen if the start and end point are the same value that the
    % first sorted value is nan. In this case take the second best.
    if isnan(sortedValues(1))
        winner = ordD(2);
    else
        winner = ordD(1);
    end
    
    % Compare the unit vector from the current end to the winning start
    % with all other unit vectors from the end to the other starts. If
    % there is another start such that the dot product between the winning
    % start and the other start is larger than some threshold
    % SAME_DIRECTION_THRESHOLD and that other start is closer to the
    % current end point, then choose the other start. See bug-report image
    % rebuildBoundary_bugFix_20161021_1637 for an example of the problem
    % this code fixes.
    
    uvWin = r(winner,:);
    possibleWinners = find((d < d(winner)) & (sum(uvWin .* r,2) > SAME_DIRECTION_THRESHOLD));
    
    if ~isempty(possibleWinners)
        [~,newwinner] = min(d(possibleWinners));
        winner = possibleWinners(newwinner);
    end
    
    
    % If the start of the current group is closer than the start of any
    % other group, then take the next closest group.
%     if winner == currentGroup
%         winner = ordD(2);
%     end
    
    % Place the winning group as the next group in the boundary.
    newOrder(i+1) = winner;
    currentGroup = winner;
end
% Save the newOrder.
if debug
    Info.newOrder = newOrder;
end

if newOrder(1) ~= newOrder(end)
    % If the first and last group are not the same, then remove both of
    % them. We remove both of them because we know the last group should be
    % a dubplicate, and we remove the first one because if we did not come
    % back to it then it is a far away group.
    newOrder([1,end]) = [];
else
    % If the first and last group are the same, then just remove the last
    % group.
    newOrder(end) = [];
end

% Save the newOrder.
if debug
    Info.finalOrder = newOrder;
end


% If the newOrder and the originalOrder are not the same, the create the
% new boundary.
if ~isequal(newOrder,originalOrder)
    
    partsb = cell(numel(groups)+1,1);
    partsn = cell(numel(groups)+1,1);
    
    for i = 1:numel(groups)+1
        partsb{i} = B(bnds(i)+1:bnds(i+1),:);
        partsn{i} = n(bnds(i)+1:bnds(i+1),:);
    end
    B = cat(1,partsb{newOrder});
    n = cat(1,partsn{newOrder});
    
    % Get the square distance between adjacent boundary vertices.
    d = sum(diff(B,[],1).^2,2);

    % Find the groups of edges that make up the polygon by looking for any two
    % adjacet boundary vertices with a distance larger than 2.
    groups = find(d > 2);

    % If there are not groups, then our boundary only had one group. We just
    % need to see if the start point and end point are also right next to each
    % other, and, if not, then add in one cut from the start to the end.
    if isempty(groups)
        if sum((B(end,:) - B(1,:)).^2) > 2
            [cuts,dcuts] = orderCutVertices(B(1,:),B(end,:),n(1,:),n(end,:));
        end
        return;
    end
end

% Initialize the cuts array.
if sum((B(end,:) - B(1,:)).^2) > 2
    cuts = zeros(numel(groups)+1,4);
    dcuts = zeros(numel(groups)+1,4);
    [cuts(end,:),dcuts(end,:)] = orderCutVertices(B(1,:),B(end,:),n(1,:),n(end,:));
else
    cuts = zeros(numel(groups),4);
    dcuts = zeros(numel(groups),4);
end

% Calculate the cuts for each group.
for i = 1:numel(groups)
    [cuts(i,:), dcuts(i,:)] = orderCutVertices(B(groups(i),:), B(groups(i)+1,:), n(groups(i),:), n(groups(i)+1,:));
end

end
% rebuildBoundary : changeLog
% 2016-10-18 : Re-wrote alignCuts()
% 2016-10-sometimes : moved alignCuts to its own function and renamed
%                     orderCutVertices()
% 2016-10-21 : Fixed bug in determining the next start. example of bug may
%              be see in bug-report image
%              rebuildBoundary_bugFix_20161021_1637.png