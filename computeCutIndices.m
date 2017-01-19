function idx = computeCutIndices(cut,cutNorm,B,n)
% COMPUTECUTINDICES  Compute the index of the cut vertices in the boundary
% B.
%
% idx = computeCutIndices(cut,cutNorm,B,n)
%
% Input parameters:
% cut - Nx4 array where cut(i,1:2) and cut(i,3:4) give the x-y locations of
%       the two vertices of the i'th cut, respectively.
% cutNorm - Nx4 array where cutNorm(i,1:2) and cutNorm(i,3:4) give the unit
%           normals (point inward) for the vertices cut(i,1:2) and
%           cut(i,3:4), respectively.
% B - Mx2 array with the boundary vertices. Holes should be delimited by a
%     row of NaNs.
% n - Mx2 array with the unit normals (pointing inward) at each boundary
%     vertex.
%
% Output parameters:
% idx - Nx2 array giving the index of of each cut vertex in the boundary B.
%       idx(i,1) and idx(i,2) give the index for vertex cut(i,1:2) and
%       cut(i,3:4), respectively.
%
% See also COMPUTEBOUNDARYINDICES

% JamesKapaldo
% 2016-10-15

cutSize = size(cut);

% Reshape the cuts to 2*N x 2
cut = reshape(cut',2,cutSize(1)*2)';
cutNorm = reshape(cutNorm',2,cutSize(1)*2)';

% Compute the indices
idx = computeBoundaryIndices(cut,cutNorm,B,n);

% Reshpae the indices to Nx2
idx = reshape(idx',2,cutSize(1))';