%% Looking at how delta beta scales when you increase the number of segments
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
% H = H(1:10, 1:10);
H = UpdateH(H, 1, 0);
nCarlo = 100;
beta0 = 2;
beta = beta0;
segmentsVec = 2.^([4:8]);

figure(30);clf;shg

final_state = [];
for ii = 1:length(segmentsVec)
    nb_segments = segmentsVec(ii);
    % scale coupling rates so total time scale matches that of 16 segments
    scale = (nb_segments/16);
    betaScale = sqrt(2)^(log2(nb_segments/16));
    tau_eff = 1/scale;
    
    diss_ops = (beta0/betaScale)^2/12 / tau_eff
    diss_ops = (beta0^2)/12

    % This is needed becase in the sim beta (is really)= beta*time
    
    
    Rho = SolveWG_Sim(H/scale, beta0/betaScale,...
        nCarlo, nb_segments);
%     beta = beta0/(nb_segments/16);
%     (nb_segments/16)
    subplot(2, 4, ii)
    pcolor(Rho)
    if ii == 1
        title(['sqrt(2) scaling (beta0:', num2str(beta0)])
    else
        title('sqrt(2) scaling')
    end
    
    shading interp
    xlabel('segs')
    ylabel('site')
    
    final_state = [final_state, Rho(:,end-3:end)];
end

figure(31);clf;shg
pcolor(final_state)
title('convergence?')
xlabel('last 3 times per sim (stacked segments)')
ylabel('site')
shading interp
%% Unitary simulation compared with analytic solution
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = H;
segments = 16;
scale = segments/16;
analytic = [];
for tt = linspace(1, 16, segments)
    analytic = [analytic, SolveME(H, 0, tt, false, 6)];
end
simulated = SolveWG_Sim(H/scale, 0, 1, segments);


figure(32);clf;shg
subplot(3, 1, 1)
pcolor(analytic)
shading interp
title('analytic zero noise solution')
xlabel('time')
ylabel('site')
colorbar 

subplot(3, 1, 2)
pcolor(simulated)
shading interp
title('simulated zero noise solution')
xlabel('time')
ylabel('site')
colorbar 

subplot(3, 1, 3)
error = abs(simulated - analytic);
% error(:,1:(segments/4)) = 0;
pcolor(error)
shading interp
title('|analytic - simulated| ')
xlabel('time')
ylabel('site')
colorbar


%% DeltaBeta simulation compared with analytic solution
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
time_scaling = 1;
H = H*time_scaling;
segments = 16;
scale = segments/16;
nCarlo = 100;
analytic = [];


%%% converting between diss ops and delta beta
deltaBeta = .3;
deltaBeta = deltaBeta*time_scaling;
dissRate = deltaBeta^2/12;
betaScale = sqrt(2)^(log2(segments/16));

%%%
for tt = linspace(1/scale, 16, segments)
    analytic = [analytic, SolveME(H, dissRate, tt, true, 6)];
    tt
end
simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);

% simulated(:,1:(segments/8)) = 0;
% analytic(:,1:(segments/8)) = 0;



figure(33);clf;shg
subplot(3, 1, 1)
pcolor(analytic)
shading interp
title('analytic dissOps solution')
xlabel('time')
ylabel('site')
colorbar 

subplot(3, 1, 2)
pcolor(simulated)
shading interp
title('simulated random phasessolution')
xlabel('segment')
ylabel('site')
colorbar 

subplot(3, 1, 3)
error = abs(simulated - analytic);
error(:,1:(segments/2)) = 0;
plot(sum(error.^2))
shading interp
title('sum((analytic - simulated)^2)')
xlabel('time')
ylabel('site')
colorbar


%% Time dependence for fixed num_segs
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
time_scaling = 1;
H = H*time_scaling;
segments = 16;
scale = segments/16;
nCarlo = 20;

betaVec = linspace(0, 10, 4);
tVec = linspace(0, 1, 16);

%%% converting between diss ops and delta beta
figure(35);clf;shg;hold on
cmap = parula(length(betaVec));

for ii = 1:length(betaVec)
    deltaBeta = betaVec(ii)
    deltaBeta = deltaBeta*time_scaling;
    dissRate = deltaBeta^2/12
    betaScale = sqrt(2)^(log2(segments/16));

    %%%
    analytic = [];
    for tt = tVec
        analytic = [analytic, SolveME(H, dissRate, tt, true, 6)];
        tt
    end
    prob = sum(analytic(8:end, :));
    plot(tVec, prob, 'color', cmap(ii, :))
    pause(0.001)
end
% simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);

xlabel('time')
ylabel('transport prob.')
title('transport prob over time at different dephasings')
legendCell = cellstr(num2str(betaVec', 'Dbeta=%-d'));
legend(legendCell)

