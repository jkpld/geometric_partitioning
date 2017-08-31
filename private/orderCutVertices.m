function [cut,cutN] = orderCutVertices(v1,v2,n1,n2)
% ORDERCUTVERTICES  Order the vertices of a cut.
%
% The vertices will be aligned such that if the cut is more horizontal, the
% vertex with the smallest x value will be first, or if the cut
% is more vertical, the vertex with the smallest y value will be first.

% James Kapaldo
% 2016-10-26

dv = abs(v1-v2);

if dv(1) <= dv(2)
    % smallest x first
    if v1(1) <= v2(1)
        cut = [v1,v2];
        cutN = [n1,n2];
    else
        cut = [v2,v1];
        cutN = [n2,n1];
    end
else
    % smallest y first
    if v1(2) <= v2(2)
        cut = [v1,v2];
        cutN = [n1,n2];
    else
        cut = [v2,v1];
        cutN = [n2,n1];
    end
end

end
% orderCutVertices : changeLog
% 2016-10-18 : Re-write, previous version commented out below
% 2016-10-26 : Re-named from createCuts to orderCutVertices


% % createCut puts the two vertices in a defined order based on their
% % distances from the origin. This is important because later on we will
% % compair cuts of differnt boundary partitions by looking at the distance
% % between the cut endpoints, and, if they were not in a defined order, the
% % calculation would be more difficult.
% % n1 and n2 are the normals at the vertex's locations, we will just put
% % them in the same order as v1,v2.
%     d1 = sum(v1.^2);
%     d2 = sum(v2.^2);
%     if (d1 == d2 && v1(1) < v2(1)) || d1 < d2
%         cut = [v1,v2];
%         cutN = [n1,n2];
%     else
%         cut = [v2,v1];
%         cutN = [n2,n1];
%     end