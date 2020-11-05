function out = GeneralHamiltonian(E0, DeltaE, C0, DeltaC, dims)

rng('shuffle');

energies = E0 + (rand(1, dims)-0.5)*DeltaE;

couplings = C0 + (rand(1, dims-1)-0.5)*DeltaC;

H0 = diag(energies);
H0 = H0 + diag(couplings, 1) + diag(couplings, -1);

out = H0;


end