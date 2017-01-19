function dVint = interparticleForce(D)

A = 1;
sig2 = 4^2;
mu = 2;

D2 = D.*D;

dVint = -1./D2 + (A*(D-mu)/sig2) .* exp(-(D2 - 2*mu*D + mu*mu)/(2*sig2));
