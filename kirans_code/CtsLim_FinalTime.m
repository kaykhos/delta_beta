%% DeltaBeta simulation compared with analytic solution
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = H;
segments = 128;
scale = segments/16;
betaScale = sqrt(2)^(log2(segments/16));
nCarlo = 200;


%%% converting between diss ops and delta beta
deltaBeta = .5;
dissRate = 1*deltaBeta^2/12;
%%%

analytic = SolveME(H, dissRate, 16, true, 6);
simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale,...
    nCarlo, segments);
simulated = simulated(:,end);

 
figure(40);clf;shg;hold on
plot(simulated, 'r')
plot(analytic, 'b')
legend('sim', 'exact')


%% Loop over delta beta and segments
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = UpdateH(H, 1);
time_scaling = 1;
H = H*time_scaling;
figure(41);clf;shg
nCarlo = 4;
segVec = 2.^[4];
betaVec = linspace(0, .5, 5);
betaVec = betaVec*time_scaling;

for ii = 1:length(segVec)
    for jj = 1:length(betaVec)
        segments = segVec(ii);
        scale = segments/16;
        betaScale = sqrt(2)^(ii-1);
        deltaBeta = betaVec(jj);
        dissRate = 1*deltaBeta^2/12;
        
        analytic = SolveME(H, dissRate, 16, true, 6);
        simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);
        simulated = simulated(:,end);
              
        subplot(length(segVec), length(betaVec), jj + (ii-1)*(length(betaVec)))
        hold on
        plot(simulated, 'r')
        plot(analytic, 'b')
        legend('wg sim', 'exact')
        title(['beta = ', num2str(deltaBeta)])
        ylabel(['segs = ', num2str(segments)],'fontweight','bold')
    end
end


%% Loop over dissRate
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = UpdateH(H, 3, 10);
time_scaling = 20;
H = H*time_scaling;
figure(41);clf;shg
nCarlo = 50;
segVec = 2.^[4:8];
betaVec = linspace(.5, .5, 1);
betaVec = betaVec*time_scaling;

for ii = 1:length(segVec)
    for jj = 1:length(betaVec)
        segments = segVec(ii);
        scale = segments/16;
        betaScale = sqrt(2)^(log2(segments/16));
        deltaBeta = betaVec(jj);
        
        dissRate = 1*deltaBeta^2/12
        dissRate2 = (deltaBeta/betaScale)^2*(1/scale)/12
        
        analytic = SolveME(H, dissRate, 16, true, 6);
        simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);
        simulated = simulated(:,end);
              
        subplot(length(segVec), length(betaVec), jj + (ii-1)*(length(betaVec)))
        hold on
        plot(simulated, 'r')
        plot(analytic, 'b')
        legend('wg sim', 'exact')
        title(['beta = ', num2str(deltaBeta/time_scaling)])
        ylabel(['segs = ', num2str(segments)])
    end
end


%% Plotting transiiton probability as a function of delta beta
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = UpdateH(H, 3, 2);
time_scaling = 1;
H = H*time_scaling;
nCarlo = 1000;
segVec = 2.^[4]; 
betaVec = linspace(0, .5, 10);
betaVec = betaVec*time_scaling;
    

transition_probability_exact = zeros(length(segVec), length(betaVec));
transition_probability_sim = zeros(length(segVec), length(betaVec));
for ii = 1:length(segVec)
    for jj = 1:length(betaVec)
        segments = segVec(ii);
        scale = segments/16;
        betaScale = sqrt(2)^(log2(segments/16));
        deltaBeta = betaVec(jj);
        dissRate = 1*deltaBeta^2/12;
        
        analytic = SolveME(H, dissRate, 16, true, 6);
        simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);
        simulated = simulated(:,end);
        
        prob_exact = sum(analytic(8:end));
        prob_sim = sum(simulated(8:end));
        
        transition_probability_exact(ii, jj) = prob_exact
        transition_probability_sim(ii, jj) = prob_sim
    end
end

figure(45);clf;shg;hold on
plot(betaVec/time_scaling, transition_probability_exact', 'r.')
plot(betaVec/time_scaling, transition_probability_sim', 'k.')
xlabel('delta beta')
ylabel('transition probability')
title('WG sim probabilieits')
legend('exact', 'sim')
