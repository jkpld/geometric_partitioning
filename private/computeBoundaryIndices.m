function idx = computeBoundaryIndices(vertex,vertexNorm,B,n)
% COMPUTEBOUNDARYINDICES Find the index of a vertex in a boundary.
%
% idx = computeBoundaryIndices(vertex,vertexNorm,B,n) will compute the
% index of the point with an x-y location given by vertex and
% unit normal (pointing inward) given by vertexNorm on the boundary B with
% normals n (inward pointing).
%
% Input parameters:
% vertex - Nx2 array of x-y locations on the boundary B
% vertexNorm - Nx2 array of the unit normals (pointing inward) at each of
%              the locations in vertex
% B - Mx2 array with the boundary vertices. Holes should be delimited by a
%     row of NaNs.
% n - Mx2 array with the unit normals (pointing inward) at each boundary
%     vertex.
%
% Output parameters:
% idx - Nx1 array giving the index of vertex in B.
%
% Notes: The boundary B should be the boundaries of an object from a BW
% mask. This means that the boundary cannot cross itself, but it can pass
% through the same x-y points. Where it does pass through the the same x-y
% points the normals will be pointed in oposite directions.
%
% See also COMPUTEBOUNDARYINFORMATION

% JamesKapaldo
% 2016-10-15

% Offset the boundary by the normal vector so that there are no repeat
% locations in the boundary
offset_B = B + 0.25*n;

% Offset the vertex by the unit normals so that they overlap with the
% offset boundary
offset_vertex = vertex + 0.25*vertexNorm;

% Find the indices of the cut vertices.
[Bidx,vrtIdx] =  find( (offset_B(:,1) == offset_vertex(:,1)') & (offset_B(:,2) == offset_vertex(:,2)') );

% Correct for vertices that are not on the boundary (in the case of
% triangle centers)

idx = nan(size(vertex,1),1);
idx(vrtIdx) = Bidx;

end
% computeBoundaryIndices : changeLog
% 2016-10-23 : fixed bug from NaN vertex/vertexNorm, which are from
%              triangle centers.