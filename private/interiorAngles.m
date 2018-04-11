function th = interiorAngles(B)
% INTERIORANGLES  Calculate the interior angles of a polygon given by B
% where each row is a vertex.
%
% The polygon must be in clockwise orientation, otherwise the exterior
% angles will be returned.
%
% The interior angle is compute by atan2(crossproduct,dotproduct).

% James Kapaldo
% 2016-10-17

% http://stackoverflow.com/questions/28821329/interior-angles-of-irregular-polygon-with-angles-180

if ~isPolyCW(B)
    B = flipud(B);
end

v1 = circshift(B,1,1) - B;
v2 = circshift(B,-1,1) - B;

th = atan2(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1), sum(v1.*v2,2));
th(th<0) = th(th<0) + 2*pi;

end

function tf = isPolyCW(xy)

if sum((xy(1,:) - xy(end,:)).^2) < 1e-3
    xy(end,:) = [];
end

dups = [false; any(diff(xy,1,1)==0,2)];
xy(dups,:) = [];

xy(:,1) = xy(:,1) - mean(xy(:,1));
n = size(xy,1);
if n < 3
    tf = true;
    return;
else
    a = sum( xy([2:n,1],1) .* (xy([3:n,1,2],2) - xy(:,2)) );
    tf = (a <= 0);
end

end
