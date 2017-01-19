function [E, dE] = electronSystemEnergy_onlyPos(y,V,dVx,dVy,Vint_fun,dVint_fun,pdistInds)
% y(1) x position of electron 1
% y(2) y position of electron 1
%
% y(3) x position of electron 2
% y(4) y position of electron 2
%   .
%   .
%   .
%   .
%
% The number of electrons will be given by N

% Number of inputs (4 times the number of particles)
N = numel(y)/2;

r = reshape(y,2,N)';

% Distance between particles
D = pdist(r);
Vint = Vint_fun(D);
% dVint = dVint_fun(D);


E = sum(V(r)) + sum(Vint);

% Direction of force.
% r_hat = r(pdistInds(:,1),:) - r(pdistInds(:,2),:);
% r_hat = r_hat ./ sqrt(sum(r_hat.^2,2));
% 
% r
% r_hat
% dVint
% Vint
% 
% % Add in all of the interaction terms dv
% dE = zeros(N,2);
% for i = 1:size(pdistInds,1)
%     dE(pdistInds(i,:),:) = dE(pdistInds(i,:),:) - dVint(i)*(r_hat(i,:).*[1;-1]);
% end
% dE
% % Add in potential term to dv
% dE = dE + [dVx(r), dVy(r)];
% dE
% dE = dE';
% dE = dE(:);