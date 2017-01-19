% PARTITIONFUNCTIONS
%
% Files
%   calculateCuts                        - calculateCuts will compute the cuts to boundary B when each boundary
%   checkCutIntegrity                    - Determine if the cuts intersect the boundary and
%   chooseWinningCuts                    - Decide between the triangle cuts and the cuts
%   computeBoundaryIndices               - Find the index of a vertex in a boundary.
%   computeCutAdjacency                  - Determine which cuts are adjacent to each
%   computeCutIndices                    - Compute the index of the cut vertices in the boundary
%   constrainedObjectCenterTriangulation - Compute constrained triangulation
%   createObjectImages                   - Create the mask, edge, and intensity image for the
%   createObjectPartitions               - Partition a boundary into disjoint regions.
%   createTriangleCuts                   - Create cuts from the from the centers of the
%   findAssociatedCuts                   - Determine the cuts associated with each triangle
%   getSearchInds                        - return the indices within searchRange of ind. 
%   interiorAngles                       - Calculate the interior angles of a polygon given by B
%   optimizeCuts                         - Search in the local region of the vertices that make up a
%   orderCutVertices                     - Order the vertices of a cut.
%   rebuildBoundary                      - Take in the pieces of boundary assigned to a center and
%   triGroups                            - This will recursively find all groups of triangles that are joind by an
