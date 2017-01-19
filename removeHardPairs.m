function toRemove = removeHardPairs(cut,cutK,cutAdjacency)
% REMOVEHARDPAIRS  Remove one of the cuts in a pair of cuts that are
% adjacent at both vertices. Keep the cut with the largest curvature and
% shortest distance.

% James Kapaldo
% 2016-10-26

toRemove = false(size(cut,1),1);

A = triu(cutAdjacency);
[i,j] = find(A==2);

for pairNum = 1:numel(i)
    s1 = sum(cutK(i(pairNum),:)) / sqrt(sum((cut(i(pairNum),1:2)-cut(i(pairNum),3:4)).^2,2));
    s2 = sum(cutK(j(pairNum),:)) / sqrt(sum((cut(j(pairNum),1:2)-cut(j(pairNum),3:4)).^2,2));

    if s2 < s1
        toRemove(j(pairNum)) = true;
    else
        toRemove(i(pairNum)) = true;
    end
end




end