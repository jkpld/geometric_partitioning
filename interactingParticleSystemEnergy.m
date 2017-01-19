%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                VECTORIZED VERSION OF CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function E = interactingParticleSystemEnergy(t,y,extraInputs)
% t is the time
% y should accept a matrix of inputs where each row is a set of
% variables to calculate the derivatives of
%
% y = [y1 | y2 | y3 | ... | yn]';
%
% y1 = [ rx1, ry1, vx1, vy1, rx2, ry2, vx2, vy2, ... , rxN, ryN, vxN, vyN ]'
%
%
y = y';
% Get all the inputs needed ==============================================
% Particle properties. 
% These could be one of the following types
% * An array of Nx1
% * A function handle that takes in one parameter, the time.
% * A griddedInterplant that takes in one parameter, the time.
q = extraInputs.q; % particle charges (This charge should also include the square root of the coupling constant -- k or 1/(4*pi*eps0).)
m = extraInputs.m; % particle masses 
% alpha = extraInputs.alpha; % particle damping coefficients 

if isa(q,'griddedInterpolant') || isa(q,'function_handle')
    q = q(t);
end
if isa(m,'griddedInterpolant') || isa(m,'function_handle')
    m = m(t);
end
size(m)
% if isa(alpha,'griddedInterpolant') || isa(alpha,'function_handle')
%     alpha = alpha(t);
% end

% Field properties and particle interactions.
V = extraInputs.V; % griddedInterpolant giving V
Vint_fun = extraInputs.Vint; % function handle for particle intercation - takes in distance between two particles an outputs a scalar

% Indices for particle pairs.
pdistInds = extraInputs.pdistInds; % N*(N-1)/2 x 2 - indices linking each distance between two partices returned by pdist to the two particles.
accumIndsOut = extraInputs.accumIndsOut;
accumIndsIn = extraInputs.accumIndsIn;

% Number of inputs
[N,M] = size(y); 
% N = 4 * (number of particles)
% M = number of different input vectors

% Offset indices
offset = (0:4:N-1);

N = N/4;

% Indices of position and velocity for each particle
rInds = [1;2] + offset;
pInds = [3;4] + offset;

r = y(rInds(:),:);
p = y(pInds(:),:);

% r is now a matrix that looks like
%  [ rx1_1, rx1_2, ..., rx1_M ;
%  [ ry1_1, ry1_2, ..., ry1_M ;
%  [ rx2_1, rx2_2, ..., rx2_M ;
%  [ ry2_1, ry2_2, ..., ry2_M ;
%  [   .    .             .
%  [   .       .          .
%  [   .          .       .
%  [ rxN_1,      ...    rxN_M ]
%  [ ryN_1,      ...    ryN_M ]
%
% and v looks the same

% In order to calculate the potential at each nuclei position, we need the
% position vectors to be in the form 
%
% [ x1,y1; 
%   x2,y2; 
%   x3,y3; 
%    ... ; 
%   xN*M,yN*M]
%
% (This is the form expected by the griddedInterpolant function.)

r = reshape(r,[2,N*M])';
E = V(r); % N*M x 1

% In order to calculate the interaction between particles, we will reshape
% each set of particles to a page. Thus it will have N rows with two
% columns (x, y) and M pages
% 
% Page 1
% [ rx1_1 ry1_1 ;
% [ rx2_1 ry2_1 ;
% [   .
% [   .
% [ rxN_1 ryN_1 ];
%
%  .
%  .
%  .
% 
% Page M
% [ rx1_M ry1_M ;
% [ rx2_M ry2_M ;
% [   .
% [   .
% [ rxN_M ryN_M ];

r = permute( reshape( r', [2, N, M]), [2,1,3]); % Nx2xM

% Convert dp into the same form.
E = permute( reshape( E', [1, N, M]), [2,1,3]); % Nx1xM

% Also, put p into the same form
p = permute( reshape( p, [2, N, M] ), [2,1,3]); % Nx2xM

% Reshape q and m
if numel(m)~=1
    if numel(m)~=N
        m = reshape(m,N,1,M);
    end
end
if numel(q)~=1
    if numel(q)~=N
        q = reshape(q,N,1,M);
    end
end

% Get the distance between particles and save them as an N*(N-1)/2x1xM
% matrix
D = zeros(N*(N-1)/2,1,M);
for i = 1:M
    D(:,1,i) = pdist(r(:,:,i));
end

% Get the interaction force between the particles (which only depends on
% the distance between them)
Vint = Vint_fun(D); % N*(N-1)/2 x 1 x M

for i = 1:M
    E(:,:,i) = E(:,:,i) + (accumIndsOut*Vint(:,:,i) + accumIndsIn*Vint(:,:,i))/2;
end


% Add the interaction terms into dp
% for i = 1:size(pdistInds,1)
%     E(pdistInds(i,:),:,:) = E(pdistInds(i,:),:,:) + prod(q(pdistInds(i,:))) * Vint(i,1,:);
%     % Notes
%     %
%     % pdistInds is N*(N-1)/2 by 2 matrix. (returned by getPdistInds(N))
%     %
%     % E(pdistInds(i,:),:,:) will be 2x1xM. Each row will be the two
%     % particles interacting with each. 
%     %
%     % Vint(i,1,:) is a 1x1xM.
%     %
%     % This for loop could probably be removed by some reworking of the
%     % arrays and a call to accumarray, but I'm not sure how at the moment
%     % (2016-09-28).
% end

% Do not consider damping
% E = E - (alpha./m) .* p;

% Add in the kinetic energy
E = E + sum(p.^2,2) ./ (2*m); % Nx1xM

% Reshape E to same form as input (M sets by N particles)
E = squeeze(E)'; % MxN







