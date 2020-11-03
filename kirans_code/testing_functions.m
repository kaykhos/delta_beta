%% Testing matrix vs vectorised versions of the unitary solution

Generated8SiteHL_sqrtC_decoherence_and_dephasing;

time = 16;
init = 6;
dissRate = 0;

figure(101);clf;shg
for init = 1:6
    mat1 = SolveME(H, dissRate, time, false, init);
    vec1 = SolveME(H, dissRate, time, true, init);
    
    subplot(2, 3, init);hold on
    plot(mat1)
    plot(vec1, 'd')
    legend('mat', 'vec')
    xlabel('site')
    ylabel('prob')
    title(['init state ', num2str(init)])
end


%% movie of dissipative hopping
figure(102);clf;shg

for tt = 1:3:40
    vec1 = SolveME(H, 1, tt, true, 6);
    mat1 = SolveME(H, 0, tt, false, 6);

    subplot(2, 1, 1)
    plot(vec1);
    title(['t = ', num2str(tt), ': diss'])
    axis([1, length(vec1), 0, 0.2])
    
    subplot(2, 1, 2)
    plot(mat1);
    title(['t = ', num2str(tt), ': unitary'])
    axis([1, length(vec1), 0, 0.2])
    shg
    pause(0.001)
end





%% Loop over dissRate
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
% H = H(1:8,1:8);
figure(103);clf;shg
nCarlo = 30;
segVec = 2.^[4:12];
betaVec = linspace(.01, .33, 1);

cmap = parula(length(segVec));
for jj = 1:length(betaVec)
    convergence_mat = [];
    for ii = 1:length(segVec)
        segments = segVec(ii);
        scale = segments/16;
        betaScale = sqrt(2)^(ii-1);
        deltaBeta = betaVec(jj);
        
        simulated = SolveWG_Sim(H/scale, deltaBeta/betaScale, nCarlo, segments);
        simulated = simulated(:,end);
        convergence_mat = [convergence_mat, simulated]; 
        
        subplot(1, length(betaVec), jj);hold on
        plot(simulated, 'color', cmap(ii,:))

    end
%     subplot(1, length(betaVec), jj)
%     pcolor(convergence_mat)
%     colorbar
    shading interp
    title(['beta = ', num2str(deltaBeta)])
    ylabel('prob')
    xlabel('site')
    legend(cellstr(num2str(segVec', 'N=%-d')));
end



%% What does delta beta do to the simulations? 
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = UpdateH(H, 3, 0);
time_scale_vec = linspace(0, 50, 16);

nCarlo = 20;
nb_segments = 256;
scale = nb_segments/16;


beta0 = 1;
betaScale = sqrt(2)^(log2(nb_segments/16));


diss_sim = [];
unit_sim = [];
for jj = 1:length(time_scale_vec)
    tt = time_scale_vec(jj);


    simulated = SolveWG_Sim(H*tt/scale, beta0*sqrt(tt)/betaScale,...
        nCarlo, nb_segments);
    
    unitary = SolveWG_Sim(H*tt/scale, 0, nCarlo, nb_segments);
    
    diss_sim = [diss_sim, simulated(:,end)];
    unit_sim = [unit_sim, unitary(:,end)];
end


figure(105);clf;shg
for ii = 1:length(time_scale_vec)
    subplot(4, 4, ii);hold on 
    plot(diss_sim(:,ii), '.')
    plot(unit_sim(:,ii))
    legend('diss', 'unit')
    title(num2str(time_scale_vec(ii)))
end
