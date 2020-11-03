function out = SolveWG_Sim(H, deltaBeta, nCarlo, segments, zero_mean)
%NewSimulation new simulatoin of a waveguide set
%   Based on Hao's code but modified for speed

n=size(H,1);
out=zeros(n,segments);  %

for calc=1:nCarlo
    Psi=zeros(n,1);  %initial state for each MC simulation
    Psi(6)=1;% 
    if zero_mean
        R=deltaBeta.*(rand([segments,7]) - 0.5); % random number generation
        % rows = each segment of WG, only first 7 guids are dynamic
    else
        R=deltaBeta.*(rand([segments,7])); % random number generation
    end


    for ii=1:1:segments
        H_beta = diag([R(ii,:), zeros(1, n-7)]); % construch deltaBeta H
        H_t = H + H_beta; % Total Hamiltonian for this seg
        U = expm(-1i*H_t); % Unitary for 1mm Only)
        Psi = U*Psi; % update wave function
        Prob_t = abs(Psi.^2); % calculate and store probability
        out(:,ii) = out(:,ii) + Prob_t;
    end
end
out = out / nCarlo;


% print statistics of random phase (useful for looking at convergence of
% nCarlo)
phases = deltaBeta*(rand([segments, nCarlo])-.5);    
total_phase = sum(phases);
disp(['std of total accumilated phase = ', num2str(std(total_phase))])

end