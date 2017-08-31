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

if ~ispolycw(B(:,1),B(:,2))
    B = flipud(B);
end

v1 = circshift(B,1,1) - B;
v2 = circshift(B,-1,1) - B;

th = atan2(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1), sum(v1.*v2,2));
th(th<0) = th(th<0) + 2*pi;

end