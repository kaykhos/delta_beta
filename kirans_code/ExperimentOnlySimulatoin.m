%% Actual time dependence (full ME check)
H = GenHWithUnits(30, 1.3);
H = UpdateH(H, .5, 10); % try to change eigen basis

nb_segments = 16;
nCarlo = 1;
% betaVec = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
betaVec = linspace(0, 1, 14);
nb_wg_blocks = 100;


figure(60);clf;shg
experiment_data = [];
for jj = 1:nb_wg_blocks
    transport_probability = [];
    for ii = 1:length(betaVec)
        beta = betaVec(ii);
        out = SolveWG_Sim(H, beta, nCarlo, nb_segments, true);
        final_state = out(:,end);   

        if nb_wg_blocks <= 10
            subplot(2, 7, ii);hold on
            plot(final_state)
            xlabel('site')
            ylabel('prob')
            title(['single block Dbeta = ', num2str(beta)])
        end

        transport_probability = [transport_probability, ...
                                 sum(final_state(8:end))];
    end
    experiment_data = [experiment_data; transport_probability];
end


figure(61);clf;shg;hold on
cmap = jet(nb_wg_blocks);
for ii = 1:nb_wg_blocks
    plot(betaVec, experiment_data(ii,:), '.', 'color', cmap(ii,:))
end
errorbar(betaVec, mean(experiment_data), std(experiment_data), 'kd', 'markerfacecolor','k')
xlabel('delta beta mm^{-1}')
ylabel('light in sink')
title('direct simulation of experiment')
% axis([min(betaVec), max(betaVec), min(min(experiment_data)), max(max(experiment_data))])
% legend(cellstr(num2str((1:nb_wg_blocks)')))

%% ME simulation
if length(H(:,1)) > 50
    disp('warning this is a big matrix problem, reduce the nb of sinks')
    return
end

dissVec = zeros(1,length(betaVec));
for ii = 1:length(betaVec)
    beta = betaVec(ii);
    diss = beta^2 * 1*2  ; % t = 1 e.g. per segment
    dissVec(ii) = diss;
end


figure(62);clf;shg
transport_probability = [];
for ii = 1:length(dissVec)
    ii
    diss = dissVec(ii);
    state = SolveME(H, diss, nb_segments, true, 6); 

    figure(62)
    subplot(2, 7, ii);hold on
    plot(state)
    xlabel('site')
    ylabel('prob')
    title(['ME sol. with diss = ', num2str(diss)])
    
    figure(60)
    subplot(2, 7, ii);hold on
    plot(state, 'r.')
    xlabel('site')
    ylabel('prob')
    title(['ME sol. with diss = ', num2str(diss)])
    

    transport_probability = [transport_probability, ...
                             sum(state(8:end))];
                         
                         
end


figure(63);clf;shg;hold on
plot(dissVec, transport_probability, 'd')
xlabel('diss op : deltaBeta^2 tau/12 (mm^{-1})')
ylabel('light in sink')
title('ME simulation with same delta beta')

figure(61);hold on
plot(betaVec, transport_probability, '-d')
