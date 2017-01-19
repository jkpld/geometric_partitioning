function [searchInds,currentInd] = getSearchInds(idx,validSet,walls,searchRadius)
% GETSEARCHINDS return the indices within searchRange of ind. 
%
% Input parameters
% ind - the index we want indices around
% validSet - the set of indices for which we can search.
% walls - an array containing all the index of all walls (this
%         also must include ind)
% searchRange - positive scalar giving the radius of the search range.
%
% Output parameters
% searchInds - the valid indices around ind subject to the constraint that
%              indices cannot pass any index in walls (excluding ind). The
%              index values are wrapped to the range of
%              min(validSet),max(validSet).
% currentInd - the index of ind within walls.

% James Kapaldo
% 2016-10-06


% fprintf('========================\n')
% fprintf('-> start:getSearchInds\n')

N = numel(validSet);

% Make all indices go from 1:N
minVS = min(validSet)-1;
ind = idx-minVS;
v = walls-minVS;

% ind
% v

searchRange = -searchRadius:searchRadius;

if N < numel(searchRange)
    searchInds = 1:N;
else
    searchInds = ind + searchRange;
    if any(searchInds<1)
        searchInds = [N+min(searchInds)+1:N,1:max(searchInds)];
    elseif any(searchInds>N)
        searchInds = [min(searchInds):N,1:max(searchInds)-N+1];
    end
end

% searchInds(:)'

currentInd = v == ind;

searchbounds = [v(circshift(currentInd,-1)), v(circshift(currentInd,1))];

% searchbounds

if currentInd(1) || currentInd(end)
    searchInds(searchInds <= searchbounds(1) & searchInds >= searchbounds(2)) = [];
else
    searchInds(searchInds <= searchbounds(1) | searchInds >= searchbounds(2)) = [];
end
% searchInds(:)'

currentInd = find(currentInd);

searchInds = searchInds + minVS;

% fprintf('-> end:getSearchInds\n')
% fprintf('========================\n')


end
% getSearchInds : changeLog
% 2016-10-15 Updated to work not just 1:N but from N1:N2.