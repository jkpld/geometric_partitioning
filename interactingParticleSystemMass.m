%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                VECTORIZED VERSION OF CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function M = interactingParticleSystemMass(t,extraInputs)

m = extraInputs.m;

if isa(m,'griddedInterpolant') || isa(m,'function_handle')
    m = m(t);
end

N = extraInputs.N;

m = [ones(N,2),m,m]';
M = diag(m(:));