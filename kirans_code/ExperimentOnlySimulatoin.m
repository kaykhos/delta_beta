%% Actual time dependence (full ME check)
H = GenHWithUnits(30, 1.0);
H = UpdateH(H, 3, 10); % try to change eigen basis
H = H;
nb_segments = 16;
nCarlo = 1;
% betaVec = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
betaVec = linspace(0, 12, 6);
% betaVec = [0]
nb_wg_blocks = 4;
rng('shuffle')

figure(60);clf;shg
experiment_data = [];
for jj = 1:nb_wg_blocks
    transport_probability = [];
    for ii = 1:length(betaVec)
        beta = betaVec(ii);
        out = SolveWG_Sim(H, beta, nCarlo, nb_segments, true, 0);
        final_state = out(:,end);   

        if nb_wg_blocks <= 10
            subplot(2, 3, ii);hold on
            plot(final_state)
            xlabel('site')
            ylabel('prob')
            title(['single block Dbeta = ', num2str(beta)])
            axis([1, 37, 0, .5])

        end

        transport_probability = [transport_probability, ...
                                 sum(final_state(8:end))];
    end
    experiment_data = [experiment_data; transport_probability];
end


figure(61);clf;shg;hold on
cmap = jet(nb_wg_blocks);
for ii = 1:nb_wg_blocks
    plot(betaVec, experiment_data(ii,:), 'k--')%, 'color', cmap(ii,:))
end
errorbar(betaVec, mean(experiment_data), std(experiment_data), 'rd-', 'markerfacecolor','r')
xlabel('delta beta mm^{-1}')
ylabel('light in sink')
title('direct simulation of experiment')
% axis([min(betaVec), max(betaVec), min(min(experiment_data)), max(max(experiment_data))])
% legend(cellstr(num2str((1:nb_wg_blocks)')))

% exp_sol = [28.4, 35.0, 74.0, 87.1;...
%         52.5, 58.0, 94.4, 96.0;...
%         70.5, 71.4, 95.9, 98.2;...
%         79.7, 85.0, 98.3, 96.0;...
%         90.9, 96.8, 98.3, 97.6;...
%         96.1, 99.4, 98.4, 98.0;...
%         98.2, 98.9, 97.6, 97.7;...
%         99.2, 98.0, 98.0, 97.2;...
%         99.9, 98.9, 96.5, 94.5;...
%         95.7, 98.7, 97.0, 97.7;...
%         95.8, 98.1, 96.3, 94.7;...
%         94.9, 96.6, 96.3, 96.3;...
%         92.6, 95.6, 96.4, 96.9;...
%         91.7, 92.6, 96.2, 96.9]/100;
% 
% 
% plot(betaVec, exp_sol, 'linewidth', 2)
% axis([0, 1, .7, 1])
% 
% 
% subplot(1, 2, 1);
% plot(fs_anderson)
% xlabel('site')
% ylabel('probability')
% title('anderson')
% subplot(1, 2, 2);
% plot(fs_non)
% xlabel('site')
% ylabel('probability')
% title('non-anderson')

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
