%% Finding out fast forwarding scalings delta Beta
% checked beta-> beta/sqrt(2) works when doubling nb_segments
% Doubling the time (time_scaling) must double variance of the phase
% So as H->time_sc*H, beta->sqrt(time_sc)*beta
%       e.g. phase variance has to be linear in time, while keeping nb_segs
%       fixed
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
time_scaling = 2;
H = UpdateH(H, 3, 0);
H = time_scaling*H;

nCarlo = 400;
beta0 = 5*sqrt(time_scaling);
segmentsVec = 16:2:64;
segmentsVec = linspace(16, 16, 16);
segmentsVec = 2.^([4, 5, 6]);

final_state = [];
for ii = 1:length(segmentsVec)
    nb_segments = segmentsVec(ii);
    % scale coupling rates so total time scale matches that of 16 segments
    scale = (nb_segments/16);
    betaScale = sqrt(2)^(log2(nb_segments/16));
%     betaScale = 2^(log2(nb_segments/16));

    tau_eff = 1/scale;
    
    diss_ops = (beta0/betaScale)^2/12 / tau_eff;
    diss_ops = (beta0^2)/12;

    % This is needed becase in the sim beta (is really)= beta*time
    
    
    simulated = SolveWG_Sim(H/scale, beta0/betaScale,...
        nCarlo, nb_segments);
    
    
    final_state = [final_state, simulated(:,end)];
end

figure(200);clf;shg;hold on
cmap = winter(length(segmentsVec));
for ii = 1:length(segmentsVec)
    plot(final_state(:,ii), 'color', cmap(ii,:))
end
title('converginv as nb_segs increase')
xlabel('sites')
ylabel('prob')
led = cellstr(num2str(segmentsVec'));
legend(led)
axis([0, 30, 0, .15])

%% Actual time dependence (unitary check)
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
time_scale_vec = linspace(0, 2, 20);
H = UpdateH(H, 3, 0);

nCarlo = 10;
beta0 = 0.0;
segmentsVec = 16:2:64;
segmentsVec = linspace(16, 16, 16);
segmentsVec = 2.^[5];


figure(202);clf;shg
figure(203);clf;shg
for ii = 1:length(segmentsVec)
    nb_segments = segmentsVec(ii);
    scale = (nb_segments/16);
    betaScale = sqrt(2)^(log2(nb_segments/16));
    
    prob_sim = [];
    prob_exa = [];

    for ii = 1:length(time_scale_vec)
        tt = time_scale_vec(ii);
        simulated = SolveWG_Sim(H*tt/scale, beta0*sqrt(tt)/betaScale,...
            nCarlo, nb_segments);
        prob_sim = [prob_sim, simulated(:,end)];
        
        exact = SolveME(H, 0, 16*tt, true, 6);
        prob_exa = [prob_exa, exact];
        figure(203);
        subplot(5,4,ii);hold on
        plot(exact, 'r');hold on
        plot(simulated(:,end), 'b.')
    end
    
    figure(202)
    subplot(2, 1, 1)
    pcolor(prob_sim)
    subplot(2, 1, 2)
    pcolor(prob_exa)
    
    % This is needed becase in the sim beta (is really)= beta*time
    
    
end

%% Actual time dependence (full ME check)
Generated8SiteHL_sqrtC_decoherence_and_dephasing;
time_scale_vec = linspace(0, 10, 16);
H = UpdateH(H, 3, 0);

nCarlo = 4;
beta0 = .5;
segmentsVec = 16:2:64;
segmentsVec = linspace(16, 16, 16);
segmentsVec = 2.^[4,4,4,4,4,4,4,4,4,4,4,4];


figure(202);clf;shg
figure(203);clf;shg
cmap = winter(length(segmentsVec));
for jj = 1:length(segmentsVec)
    nb_segments = segmentsVec(jj);
    scale = (nb_segments/16);
    betaScale = sqrt(2)^(log2(nb_segments/16));
    
    
    prob_sim = [];
    prob_exa = [];

    for ii = 1:length(time_scale_vec)
        tt = time_scale_vec(ii);
        

        simulated = SolveWG_Sim(H*tt/scale, beta0*sqrt(tt)/betaScale,...
            nCarlo, nb_segments);
        prob_sim = [prob_sim, simulated(:,end)];
        
        
        tau_eff = 1/scale*tt;
        diss_ops = (beta0)^2/12
        
        exact = SolveME(H, diss_ops, 16*tt, true, 6);
        prob_exa = [prob_exa, exact];
        figure(203);
        subplot(4,4,ii);hold on
        plot(exact, 'r');hold on
        plot(simulated(:,end), '.', 'color', cmap(jj,:))
        axis off 
        title(num2str(tt))
        if ii == 1
            legend('excat', 'wg sim')
        end
    end
    
    figure(202)
    subplot(2, 1, 1)
    pcolor(log(prob_sim+10))
    subplot(2, 1, 2)
    pcolor(log(prob_exa+10))
    
    % This is needed becase in the sim beta (is really)= beta*time
    
    
end

%%
deltaBeta = 10*rand;
disp(' ')
disp(['deltaBeta = ', num2str(deltaBeta)])
phases0 = deltaBeta*(rand(1, 512)-.5);
disp(['var phases are    : ', num2str(var(phases0))])
disp(['and expected to be: ', num2str(deltaBeta^2/12)])
disp(' ')
disp('after doubling the time (b-> sqrt2 b)')
disp('')
deltaBeta = sqrt(2)*deltaBeta;
phases1 = deltaBeta*(rand(1, 50000)-.5);
disp(['var phases are    : ', num2str(var(phases1))])
disp(['and expecte to be : ', num2str(2*var(phases0))])

disp(' ')
disp(' ')
disp(' ')
disp(' ')


