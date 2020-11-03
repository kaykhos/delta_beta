%% Transport simulation with true Hamiltonian

onsite_offset_and_randomness = [1, 0.5];
couplings_offset_and_randomness = [0, 0.1];
dims = 27;

energies = onsite_offset_and_randomness(1) + ...
    rand(1, dims)*onsite_offset_and_randomness(2);

couplings = couplings_offset_and_randomness(1) + ...
    rand(1, dims-1)*couplings_offset_and_randomness(2);

H0 = diag(energies);
H0 = H0 + diag(couplings, 1) + diag(couplings, -1);

H = H0;
