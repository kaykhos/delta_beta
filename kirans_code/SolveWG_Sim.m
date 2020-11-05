function out = SolveWG_Sim(H, deltaBeta, nCarlo, segments, zero_mean, var_added_time)
%NewSimulation new simulatoin of a waveguide set
%   Based on Hao's code but modified for speed

n=size(H,1);
out=zeros(n,segments);  %
H0 = H;
for calc=1:nCarlo
%     H(1:7,1:7) = (1+.5*rand)*H0(1:7,1:7);
    Psi=zeros(n,1);  %initial state for each MC simulation
    Psi(6)=1;% 
    time = 1+var_added_time*rand;

    if zero_mean
        R=deltaBeta.*(rand([segments,7]) - 0.5); % random number generation
        % rows = each segment of WG, only first 7 guids are dynamic
    else
        R=deltaBeta.*(rand([segments,7])); % random number generation
    end
%     R = 1.5*R

    for ii=1:1:segments
        H_beta = diag([R(ii,:), zeros(1, n-7)]); % construch deltaBeta
        H_t = H + H_beta; % Total Hamiltonian for this seg
         
%         for gep=1:1:6
%             for geq=gep+1:1:7
%                 if H_t(gep, geq)>0.15 %&& (gep==(geq-1))
% %                     [gep,geq]
%                     H_t(gep,geq)=sqrt(H_t(gep,geq)^2+(R(ii,gep)-R(ii,geq))^2/4);
%                     H_t(geq,gep)=H_t(gep,geq);
%                 end
%             end
%         end

        U = expm(-1i*H_t*time); % Unitary for 1mm Only)
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