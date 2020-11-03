function H = UpdateH(H, onsite_randomness, seed)
%GENERATEHAMILTONIAN generates a hamiltonian that has ~6 oscillatoins per
%t=1 matching the exact Hamiltonain
%   Detailed explanation goes here

if seed ~= 0
    rng(seed);
end

dims = length(H(:,1));
energies = 2*onsite_randomness*(rand(1, 7) - 0.5);
energies = [energies, zeros(1, dims - 7)];

H0 = diag(energies);
H = H + H0;

end

