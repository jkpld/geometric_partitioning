function J = interactingParticleSystemJacobian(t,y,extraInputs)
% NOTE: This function computes the Jacobian assuming Coulombic
% interactions!
%
% NOTE: this function is not vectorized to take in more than one set of y
% values.
% fprintf('in jacobian\n')
% Get all the inputs needed ==============================================
% Particle properties. 
% These could be one of the following types
% * An array of Nx1
% * A function handle that takes in one parameter, the time.
% * A griddedInterplant that takes in one parameter, the time.
q = extraInputs.q; % particle charges (This charge should also include the square root of the coupling constant -- k or 1/(4*pi*eps0).)
m = extraInputs.m; % particle masses 
alpha = extraInputs.alpha; % particle damping coefficients 

if isa(q,'griddedInterpolant') || isa(q,'function_handle')
    q = q(t);
end
if isa(m,'griddedInterpolant') || isa(m,'function_handle')
    m = m(t);
end
if isa(alpha,'griddedInterpolant') || isa(alpha,'function_handle')
    alpha = alpha(t);
end

% Field properties and particle interactions.
dVxx = extraInputs.dVxx; % griddedInterpolant giving d2V/dxdx
dVxy = extraInputs.dVxy; % griddedInterpolant giving d2V/dxdy
dVyy = extraInputs.dVyy; % griddedInterpolant giving d2V/dydy

% Indices for particle pairs.
pdistInds = extraInputs.pdistInds; % N*(N-1)/2 x 2 - indices linking each distance between two partices returned by pdist to the two particles.

% Number of inputs (4 times the number of particles)
N = numel(y);

% Offset indices
offset = (0:4:N-1)';

% Indices of position and velocity for each particle
rInds = [1,2] + offset;
% pInds = [3,4] + offset;

% Positions and momentum of each particle
r = y(rInds);
% p = y(pInds); % note needed

% Number of particls
N = N/4;

% Distance between particles
D = pdist(r)';
D5 = D.^5;
D3 = 1./(D.^3);

% Second derivative of potential at particle locations
Vab = zeros(2,2,N);
Vab(1,1,:) = -dVxx(r);
Vab(1,2,:) = -dVxy(r);
Vab(2,2,:) = -dVyy(r);
Vab(2,1,:) = Vab(1,2,:);

% Interparticle directions.
r_diff = r(pdistInds(:,1),:) - r(pdistInds(:,2),:);

int = - 3 * (r_diff .* permute(r_diff,[1,3,2])) ./ D5;
int(:,1,1) = int(:,1,1) + D3;
int(:,2,2) = int(:,2,2) + D3;
int = prod(q(pdistInds),2) .* int; % Multiply the charge in

int = permute(int, [3,2,1]); %2 x 2 x N*(N-1)/2 : Now each page is the 2x2 interaction terms that show up in the jacobian.
% Just need to some several of the pages of int to make the diagonal block
% terms.

J = zeros(4*N,4*N);

diag_dpdr = Vab;
idx12 = [1,2];
idx34 = [3,4];
for i = 1:size(pdistInds,1)
    diag_dpdr(:,:,pdistInds(i,:)) = diag_dpdr(:,:,pdistInds(i,:)) + int(:,:,i);
    
    % diag_dpdf(:,:,pdistInds(i,:)) is 2x2x2
    % int(:,:,i) is 2x2x1
    
    % Assign the off diagonal dpdr terms
    tmp_inds_c = idx12 + (pdistInds(i,1)-1)*4;
    tmp_inds_r = idx34 + (pdistInds(i,2)-1)*4;
    
    J(tmp_inds_r,tmp_inds_c) = -int(:,:,i);
    
    tmp_inds_c = idx12 + (pdistInds(i,2)-1)*4;
    tmp_inds_r = idx34 + (pdistInds(i,1)-1)*4;
    
    J(tmp_inds_r,tmp_inds_c) = -int(:,:,i);
end

% Now assign the diagonal blocks
for i = 1:N
    J(idx12+(i-1)*4,idx34+(i-1)*4) = eye(2)/m(i);
    J(idx34+(i-1)*4,idx12+(i-1)*4) = diag_dpdr(:,:,i);
    J(idx34+(i-1)*4,idx34+(i-1)*4) = - (alpha(i)/m(i)) * eye(2);
end





