function dy = electronSystem(t,y,dVx,dVy,Fint_fun,pdistInds)
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

% Number of particls
N = N/4;

% Distance between particles
D = pdist(r)';
Fint = Fint_fun(D);

% Examples of Fin_fun
%
% Coulomb force
% F(D) = 1/D.^2 
%
% close range strong repulsion, mid range attraction, long range repulsion
% F(D) = heavyside(1-D).*(min(1./D.^3, 1e9)) - heavyside(5-D).*heavyside(D-1).*(1./D.^2) + heavyside(D-5).*(1/D.^2))

% Direction of force.
r_hat = r(pdistInds(:,1),:) - r(pdistInds(:,2),:);
r_hat = r_hat ./ sqrt(sum(r_hat.^2,2));

% Add in all of the interaction terms dv
dv = zeros(N,2);
for i = 1:size(pdistInds,1)
    dv(pdistInds(i,:),:) = dv(pdistInds(i,:),:) + Fint(i)*(r_hat(i,:).*[1;-1]);
end


% Add in potential term to dv
tau = 500;
dv = dv - [dVx(r), dVy(r)] - v*(t/3000);%exp(-t/tau);% - heaviside(t-1000)*v*exp(-(t-1000)/tau);%(5*2*pi);

% Assigned dr and assume speeds are exponentially damped.
% tau = 200;

% dr = v * exp(-t/tau);
dr = v;

% Put in proper order
dy = zeros(size(y));
dy(rInds) = dr;
dy(vInds) = dv;




