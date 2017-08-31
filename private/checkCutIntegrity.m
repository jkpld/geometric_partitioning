function boundaryIntersect = checkCutIntegrity(cuts,B,useGPU)
% CHECKCUTINTEGRITY  Determine if the cuts intersect the boundary and
% output boolian: true if they intersect, false if they do not.

% James Kapaldo
% 2016-10-28


% Get cut unit vectors
uv = cuts(:,1:2) - cuts(:,3:4);
uv = uv ./ sqrt(sum(uv.^2,2));

% Offset cuts by 0.5 along their normal vector so that the end points do
% not intersect with the boundary.
cuts = cuts + 0.5*[-uv, uv];

cutSegments = nan(3*size(cuts,1),2);
cutSegments(1:3:end) = cuts(:,1:2);
cutSegments(2:3:end) = cuts(:,3:4);

% [~,~,i1] = intersections(cutSegments(:,1),cutSegments(:,2));
% if any(i1)
%     selfIntersect = true;
% else
%     selfIntersect = false;
% end

[~,~,i1] = intersections(cutSegments(:,1),cutSegments(:,2),B(:,1),B(:,2),0,0,useGPU);
i1 = i1-floor(i1);
if any(i1)
    boundaryIntersect = true;
else
    boundaryIntersect = false;
end

end