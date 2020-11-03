%% General ideas about these scripts
% The base time is 16, and should be increased by 
%       powers of two
% To keep the time scale fixed to 16 steps use, 
% H -> H/scale (where scale ensures H*max_segs = H*16)
% DeltaBeta -> DeltaBeta


%% plotting the eigenstates of the 27 site 
for ii = 1:100

    H = GenHWithUnits(20);
    % GeneralHamiltonian
    H = UpdateH(H, 5, ii);
    [vecs, vals] = eigs(H + 3*eye(size(H)), 27);
    vals = diag(vals);

    figure(99);clf;shg
    subplot(2, 2, 1);hold on
    plot(vecs(:,1:7))
    plot([6, 6], [-0.8, 0.8], 'r')
    legendCell = cellstr(num2str(vals(1:7)));
    legend(legendCell)
    title('Eigenvectors')
    xlabel('site')
    ylabel('amplitude')

    subplot(2, 2, 2);hold on
    plot(vecs(:,8:14))
    plot([6, 6], [-0.8, 0.8], 'r')
    legendCell = cellstr(num2str(vals(8:14)));
    legend(legendCell)

    subplot(2, 2, 3);hold on
    plot(vecs(:,15:21))
    plot([6, 6], [-0.8, 0.8], 'r')
    legendCell = cellstr(num2str(vals(15:21)));
    legend(legendCell)
    ylabel('amplitude')

    subplot(2, 2, 4);hold on
    plot(vecs(:,22:end))
    plot([6, 6], [-0.8, 0.8], 'r')
    legendCell = cellstr(num2str(vals(22:end)));
    legend(legendCell)
    tmp = input('good?')
end

%% Plotting hopping dynamics as a function of DeltaBeta

Generated8SiteHL_sqrtC_decoherence_and_dephasing;
GeneralHamiltonian
Hred = H;

betaVec = linspace(0, 1, 4);
nb_segments = 16;
nCarlo = 150;
figure(1);clf;shg
for ii = 1:4
    Rho = SolveWG_Sim(Hred, betaVec(ii), nCarlo,nb_segments);
    
    subplot(2, 2, ii)
    pcolor(Rho)
    title(['\Delta \beta ', num2str(betaVec(ii))])
    xlabel('t')
    ylabel('site')
    shading interp
end





%% Increasing nb_points and decreasing delta beta
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
GeneralHamiltonian
Hred = H;

nb_steps = 2.^[4:11];
nCarlo = 30;
beta0 = 1  ;
beta = beta0;   

for ii = 1:length(nb_steps)
    steps = nb_steps(ii);
    H = Hred/(2^(ii-1)); % need to reduce hamiltonian to keep total time constant
    
    rho = SolveWG_Sim(H, beta, nCarlo, steps);
    beta = beta /  sqrt(2);
    figure(51)
    subplot(2, 4, ii)
    pcolor(rho)
    title(num2str(steps))
    xlabel('t')
    ylabel('site')
    shading flat
end
    

%% Showing how beta need to change so the variance in total phase is kept the same
figure(33);clf;shg

nb_steps = 2.^(4:13);
beta = 1;
nCarlo = 10000;

for ii = 1:length(nb_steps)
    steps = nb_steps(ii);
    beta
    phases = beta*(rand([nCarlo, steps])-.5);
    beta = beta /sqrt(2);
    
    total_phase = sum(phases');
    
    subplot(2, 2, 1);hold on
    plot(ii, mean(total_phase), 'd')
    title('average phase: ')
    ylabel('div by srqt(2)', 'fontweight', 'bold')

    subplot(2, 2, 2);hold on
    plot(ii, std(total_phase), 'd')
    title('std phase: ')
%     ylabel('div by srqt(2)')
end

beta = 1;
for ii = 1:length(nb_steps)
    steps = nb_steps(ii);
    beta
    phases = beta*(rand([nCarlo, steps])-.5);
    beta = beta / 2;
    
    total_phase = sum(phases');
    
    subplot(2, 2, 3);hold on
    plot(ii, mean(total_phase), 'd')
    xlabel('log_2(segments/16)')
    ylabel('div by 2', 'fontweight', 'bold')

    subplot(2, 2, 4);hold on
    plot(ii, std(total_phase), 'd')
    xlabel('log_2(segments/16)')
%     ylabel('div by 2')

    
end


    
    
    
    
   
    
    
    
    