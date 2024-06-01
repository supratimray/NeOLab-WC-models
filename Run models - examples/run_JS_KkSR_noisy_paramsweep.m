% Run ISN(KkSR) model with 21x21 mean inputs for 2 iterations for 2
% (thetamultiplier, sigmamultiplier) values (1/0,0) and (1,1) (Poissonlike)

basefolder = pwd();
savedir = fullfile(basefolder, 'ConstContrast_JS_KkSR_allinputCombs'); mkdir(savedir);
%%
niterations = 2;
iters = 1:niterations;

thetasigma_multipliers = [[1/0, 0]; [1,1]]; % Row 1: JS2014 (OU time-constant = inf, i.e. ; noise_variance = 0), Row 2: JS_KkSR2023
uniqmultids = 1:size(thetasigma_multipliers,1);

uniqE0 = 0:0.5:20;
uniqI0 = 0:0.5:20;

%% unwrap parameter combinations
t = -0.3:0.5e-4:1; % seconds - for noisy inputs
paramnames = {'E0', 'I0', 'multiplierIDs', 'iterIDs'};
params = { uniqE0, uniqI0, uniqmultids, iters};
[varVals, varSels, nsims] = unwrapParameters(params);

E0 = varVals(:,1);
I0 = varVals(:,2);
iterIDs = varVals(:,4);

multIDs = varVals(:,3);

thetamultipliers = {thetasigma_multipliers(multIDs,1),thetasigma_multipliers(multIDs,1)};
sigmamultipliers = {thetasigma_multipliers(multIDs,2),thetasigma_multipliers(multIDs,2)};
%% save directories
inputfilename = fullfile(savedir,'inputs.mat');
simulationfilename = fullfile(savedir,'simulationResults.mat');

%%
inputVal = {E0, I0};
baselineinputAmp = {0*E0, 0*I0};

ipdc = {0*E0, 0*I0}; % dc offset value for sinusoid

% Consider 2 intervals : [baseline, stimulus]
intervals = [-1/0, 0, 1/0];
% DC component
inputDC = {[],[]}; % {1}-E, {2}-I
% Cosine component - to be set to 0 amplitude
inputPk2Base = {[],[]}; 
inputFreq = {[],[]}; 
inputOnsetPhase = {[],[]}; 
for i = 1:numel(inputPk2Base) % {1}-E, {2}-I
    inputDC{i} = baselineinputAmp{i}*[1 0] + inputVal{i}*[0 1]; % specifying DC component of input for each input combination in each interval ([1,0]-baseline; [0,1]-stimulus)
    inputPk2Base{i} = [0,0];
    inputFreq{i} = [0,0]; % constant contrast stimuli; frequency of cosine signal input
    inputOnsetPhase{i} = [0,0]; % initial phase of cosine contrast mask at onset instant (start of interval)
end

[input_t, input] = generate_OUinput_rectifiedsine(t, intervals, nsims, inputDC, inputPk2Base, inputFreq, inputOnsetPhase, thetamultipliers, sigmamultipliers); 

save(inputfilename, 'simulationfilename','thetamultipliers', 'sigmamultipliers', 'intervals','inputDC', 'inputPk2Base','inputFreq', 'inputOnsetPhase', 'input_t','input','-v7.3');
%% model instantiation

JSpop = ISN_KkSR_JS_OUip(nsims, inputfilename);

JSpop.input([], input_t); % set timestamps of simulation

%% model simulation and lfp proxy
solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
JSpop.ode(solver);

t = JSpop.EIpairs.t;
rE = JSpop.EIpairs.R(1:end/2, :);
rI = JSpop.EIpairs.R(end/2+1:end, :);
lfp = - rE - rI;
save(simulationfilename, 'JSpop','lfp','inputfilename','niterations','uniqE0','uniqI0',"thetasigma_multipliers");
clear JSpop
%% Plots
thetaids = 1:2;
thetaids_names = {'JS2014', 'KkSR2022'};
for i_thsel=1:numel(thetaids)
    thetaidsel = thetaids(i_thsel); 
    modelname = thetaids_names{i_thsel};
    
    %% spectral analysis and plot
    [E0_grid,I0_grid] = meshgrid(uniqE0,uniqI0);
    selid = findIndexByParamComb(varSels,params,{[],[],thetaidsel,[]});
    t_sel_analysis = t>=0.25 & t<=1;
    
    eFR_trialwise = mean(rE(selid,t_sel_analysis),2);
    iFR_trialwise = mean(rI(selid,t_sel_analysis),2);

    lfp_fft_trialwise = abs(fft(lfp(selid,t_sel_analysis),[],2));
    E0_trialwise = varVals(selid,1); I0_trialwise = varVals(selid,2);
    
    lfp_fft = zeros([size(E0_grid), size(lfp_fft_trialwise,2)])/0; 
    eFR = zeros(size(E0_grid))/0;
    iFR = zeros(size(E0_grid))/0;
    E0_analysis = zeros(size(E0_grid))/0;
    I0_analysis = zeros(size(E0_grid))/0;

    for rowid = 1:size(E0_grid,1)
        for colid = 1:size(E0_grid,2)
            E0_analysis(rowid,colid) = E0_grid(rowid,colid);
            I0_analysis(rowid,colid) = I0_grid(rowid,colid);

            E0I0sel = (E0_trialwise == E0_analysis(rowid,colid)) & (I0_trialwise == I0_analysis(rowid,colid));
            if sum(E0I0sel(:))>0
                eFR(rowid,colid) = mean(eFR_trialwise(E0I0sel),'all');
                iFR(rowid,colid) = mean(iFR_trialwise(E0I0sel),'all');
                lfp_fft(rowid,colid,:) = mean(lfp_fft_trialwise( E0I0sel,:),1);
            end
        end
    end
    freqs = (0:sum(t_sel_analysis)-1)/(max(t(t_sel_analysis)) - min(t(t_sel_analysis))); % Hz
    gammaband = [30,80]; % Hz
    sel_gammabandfreqs = freqs>=gammaband(1) & freqs<=gammaband(2);
    freqs_gammaband = freqs(sel_gammabandfreqs);
    lfp_fft_gammaband = lfp_fft(:,:,sel_gammabandfreqs);

    [pkpower, pkfreq_id] = max(2*log10(lfp_fft_gammaband), [], 3);
    pkfreq = freqs_gammaband(pkfreq_id);

    % plotting
    figure_allinputs = figure('windowstate','maximized','InvertHardcopy','on');

    subplot(2,2,1); % e mean firing rates
    pcolor(E0_analysis, I0_analysis, eFR); shading interp;
    colormap jet; colorbar; clim([0 1]);
    xlim([min(E0_analysis(:)),max(E0_analysis(:))]);
    ylim([min(I0_analysis(:)),max(I0_analysis(:))]);
    ylabel('I_I');
    title('mean E firing rate');

    subplot(2,2,2); % i mean firing rates
    pcolor(E0_analysis, I0_analysis, iFR); shading interp;
    colormap jet; colorbar; clim([0 1]);
    xlim([min(E0_analysis(:)),max(E0_analysis(:))]);
    ylim([min(I0_analysis(:)),max(I0_analysis(:))]);
    title('mean I firing rate');

    subplot(2,2,3); % gamma band peak power
    pcolor(E0_analysis, I0_analysis, pkpower); shading interp;
    colormap jet; colorbar; clim([prctile(pkpower,10,'all'), prctile(pkpower,90,'all')]);
    xlim([min(E0_analysis(:)),max(E0_analysis(:))]);
    xlabel('I_E');
    ylim([min(I0_analysis(:)),max(I0_analysis(:))]);
    ylabel('I_I');
    title({'Gamma band peak Power', '(log_{10}(power))'});

    subplot(2,2,4); % gamma band peak frequency
    pcolor(E0_analysis, I0_analysis, pkfreq); shading interp;
    colormap jet; colorbar; clim(gammaband);
    xlim([min(E0_analysis(:)),max(E0_analysis(:))]);
    xlabel('I_E');
    ylim([min(I0_analysis(:)),max(I0_analysis(:))]);
    title({'Gamma band peak frequency','Hz'});

    sgtitle([modelname,' model'])

    %% Example input plot - JS2014 model
    IE_ToPlot = 10;
    II_ToPlot = 5;
    rbounds = [0 1]; % range for rE and rI
    NCvariableIDs = [1,2]; % which statevars correspond to E and I (in order), for phase diagram making
    
    selid = findIndexByParamComb(varSels,params,{IE_ToPlot,II_ToPlot,thetaidsel,[]});
    rEsel = rE(selid,:);
    rIsel = rI(selid,:);

    figure_singleInput = figure('WindowState','maximized','InvertHardcopy','on');
    subplot(1,2,1);
    % run describe dynamics on @(nsims) ISN_KkSR_JS_OUip(nsims)
    % TODO 
    describeDynamics(figure_singleInput, gca, @(nsims) ISN_KkSR_JS_OUip(nsims), [IE_ToPlot,II_ToPlot], NCvariableIDs, {rbounds, rbounds});
    hold on;
    plot(rEsel(1,:),rIsel(1,:),[[0.5,0.5,0.5]],'linestyle','--','displayname','Trajectory example 1');
    scatter(rEsel(1,1), rIsel(1,1),[],[[0.5,0.5,0.5]],'marker','x')
    scatter(rEsel(1,end), rIsel(1,end),[],[[0.5,0.5,0.5]],'marker','o','filled')
    plot(rEsel(2,:),rIsel(2,:),'k','displayname','Trajectory example 2');
    scatter(rEsel(2,1), rIsel(2,1),[],'k','marker','x')
    scatter(rEsel(2,end), rIsel(2,end),[],'k','marker','o','filled')
    title({'Phase diagram', ['E_E = ',num2str(IE_ToPlot)], ['E_I = ',num2str(II_ToPlot)]}); 
    xlabel('r_E'); ylabel('r_I');
    legend;
    % lfp proxy
    ax2=subplot(3,2,2);
    plot(t, rE(selid,:),'b'); 
    hold on;
    plot(t, rE(selid,:),'b'); 
    
    legend({'r_E','r_I'});
    ylabel('population activity');
    ylim([-0.1, 0.1]);
    
    ax3=subplot(3,2,4);
    plot(t, lfp(selid,:),'k');
    % lfp proxy
    legend({'r_E','r_I'});
    ylabel('population activity');
    ylim([-0.1, 0.1]);
    xlabel('time (s)');
    linkaxes([ax2, ax3],'x');
    xlim([-0.125, 0.5]);
    % PSD of lfp proxy 
    subplot(3,2,6);
    t_sel_fft = t>=0.25 & t<=1;
    fstep = 1/(max(t(t_sel_fft))-min(t(t_sel_fft)));
    plot((0:sum(t_sel_fft)-1)*fstep, 2*log10(abs(fft(lfp(selid,t_sel_fft))),'k'));
    ylabel('log_{10}(PSD)');
    xlabel('frequency (Hz)');
    xlim([0 150]);

    sgtitle([modelname,' model'])
    
end
%% 
