function E = electronSystemEnergy(~,y,V,Vint_fun,pdistInds)
% y(1) x position of electron 1
% y(2) y position of electron 1
% y(3) x speed of electron 1
% y(4) y speed of electron 1
%
% dy(1) x speed of electron 1
% dy(2) y speed of electron 1
% dy(3) x acceleration of electron 1
% dy(4) y acceleration of electron 1
%
% y(5) x position of electron 2
% y(6) y position of electron 2
% y(7) x speed of electron 2
% y(8) y speed of electron 2
%
% dy(5) x speed of electron 2
% dy(6) y speed of electron 2
% dy(7) x acceleration of electron 2
% dy(8) y acceleration of electron 2
%
%   .
%   .
%   .
%   .
%
% The number of electrons will be given by N

% Number of inputs (4 times the number of particles)
N = numel(y);

% Offset indices
offset = (0:4:N-1)';

% Indices of position and velocity for each particle
rInds = [1,2] + offset;
vInds = [3,4] + offset;

% Positions and velocities of each particle
r = y(rInds);
v = y(vInds);

% Distance between particles
D = pdist(r);
Vint = Vint_fun(D);

% Examples of Vin_fun
%
% Coulomb potential
% V(D) = 1/D 

v2 = sum(v.^2,2);

E = (1/2)*v2 + V(r);

% Add in all of the interaction terms E
for i = 1:size(pdistInds,1)
    E(pdistInds(i,:)) = E(pdistInds(i,:)) + Vint(i);
end
