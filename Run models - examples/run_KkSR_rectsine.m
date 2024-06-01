% Run ISN(KkSR) model with 21x21 mean inputs for 2 iterations for 2
% (thetamultiplier, sigmamultiplier) values (1/0,0) and (1,1) (Poissonlike)

basefolder = pwd();
savedir = fullfile(basefolder, 'RectSineContrast_KkSR'); mkdir(savedir);
%%
niterations = 10;
iters = 1:niterations;

%% unwrap parameter combinations
t = (-1.024e4:0.5:3.072e4-0.5)*1e-4;  % seconds - for noisy inputs
paramnames = {'iterIDs'};
params = { iters};
[varVals, varSels, nsims] = unwrapParameters(params);

iterIDs = varVals;

%% save directories
inputfilename = fullfile(savedir,'inputs.mat');
simulationfilename = fullfile(savedir,'simulationResults.mat');

%%
% min and max peak values in mean input timeseries I_E and I_I
low = {2, 2}; % {1}-E, {2}-I
high = {11.5, 5.5};

% Consider 2 intervals : [baseline, stimulus]
intervals = [-1/0, 0,1, 1.3,1.7, 2.0,2.4, 2.7,3.1, 1/0];
% DC component
inputDC = { ...
    [low{1},high{1},low{1}*ones(1,numel(intervals)-3)], ...
    [low{2},high{2},low{2}*ones(1,numel(intervals)-3)]  ...
    }; % {1}-E, {2}-I
% Cosine component - to be set to 0 amplitude
inputPk2Base = { ...
    (high{1}-low{1})*[0,0,0,1,0,1,0,1,0], ...
    (high{2}-low{2})*[0,0,0,1,0,1,0,1,0]  ...
    }; 

freqs = [0,0,0,2.5,0,5,0,10,0];
inputFreq = {freqs,freqs}; 

inputOnsetPhase = { ...
    pi/2*ones(1,numel(intervals)-1), ...
    pi/2*ones(1,numel(intervals)-1)
    }; 


[input_t, input] = generate_OUinput_rectifiedsine(t, intervals, nsims, inputDC, inputPk2Base, inputFreq, inputOnsetPhase); 

save(inputfilename, 'simulationfilename', 'intervals','inputDC', 'inputPk2Base','inputFreq', 'inputOnsetPhase', 'input_t','input','-v7.3');
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
save(simulationfilename, 'JSpop','lfp','inputfilename','niterations');
clear JSpop
%% Plots
randid = randi(niterations,1); % example trial

figure('WindowState','maximized');
ax1 = subplot(3,1,1);
plot(t, mean(input(1:end/2,:),1),'b-','DisplayName','mean i_E');
hold on;
plot(t, mean(input(randid,:),1),'b-.','DisplayName','example i_E');
hold on;
plot(t, mean(input(end/2+1:end,:),1),'r-','DisplayName','mean i_I');
plot(t, mean(input(end/2+randid,:),1),'r-.','DisplayName','example i_I');
axis tight;
ylabel('Input drive (a.u.)');

ax2 = subplot(3,1,2);
plot(t, mean(rE(1:end/2,:),1),'b-','DisplayName','mean r_E');
hold on;
plot(t, mean(rE(randid,:),1),'b-.','DisplayName','example r_E');
plot(t, mean(rI(1:end/2,:),1),'r-','DisplayName','mean r_I');
hold on;
plot(t, mean(rI(randid,:),1),'r-.','DisplayName','example r_I');
axis tight;
ylabel('Firing rate (a.u.)');


ax3 = subplot(3,1,3);
plot(t, mean(lfp(1:end/2,:),1),'k-','DisplayName','mean LFP proxy');
hold on;
plot(t, mean(lfp(randid,:),1),'k-.','DisplayName','example LFP proxy');
axis tight;
xlabel('Time (s)');
ylabel('LFP proxy (a.u.)');

linkaxes([ax1,ax2,ax3],'x');

sgtitle('JS model with OU inputs (KkSR2023) for rectified sine inputs');