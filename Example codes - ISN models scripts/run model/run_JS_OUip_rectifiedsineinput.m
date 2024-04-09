function [basefolder, filename] =  run_KkSR_JS_OUip_rectifiedsineinput( niters, overWriteFlag, imax, doMP, thetas, sigmas, filename, basefolder)

% Running 'KkSR' model for the JS model used in Krishnakumaran and Ray, Cerebral Cortex 2023.

%%
if ~exist('noisetype','var')
    noisetype = 'ou';
%     noisetype = 'oulo5';
%     noisetype = 'oubp10';
end

if ~exist('niters','var')
    niters = 50;
end
if ~exist('sigmas','var')
    sigmas = [1];
end
if ~exist('thetas','var')
%     thetas = 2*pi*[4, 16, 1/0];
    thetas = []; 
    thetaE = 2*pi*16; thetaI = 2*pi*1; 
    thetas = [thetaE thetaI];
else
    thetaE = thetas(1); thetaI = thetas(2);
end
theta_passband = [thetaE,thetaI]/2/pi;

if ~exist('basefolder','var')
    basefolder = './KkSR_JS_OUip';
    basefolder = fullfile(basefolder, 'RectSineIp')
end

if ~exist('filename','var')
    filename = ['KkSR_JS_OUip','_',datestr(now,'ddmm_HHMM')];
end

if ~exist('overWriteFlag','var')
    overWriteFlag = 1==1; % Simulate model regardless of whether simulation results file is already found and overwrite them
end
overWriteFlag1 = overWriteFlag || ~exist(fullfile([basefolder,'simdata',filename,'.mat']),'file'...
    );


if ~exist('imax','var')
    imax = [11.5,5.5];
end

iemax = imax(1);
iimax = imax(2);

ebase = 2; ibase = 2;

%% Input and Simulation time description - MODIFY AS PER REQUIREMENTS
% simulation time step configuration
tspan = (-1.024e4:0.5:3.072e4-0.5)*1e-4; 

% input description
stimwin = [-1e10, 0, 0+0.4, 0.7, 0.7+0.4, 1.4, 1.4+0.4, 2.1, 2.1+0.4, 1e10];

e0 = [ebase, iemax, ebase, iemax, ebase, iemax, ebase, iemax, ebase, ebase]; lE=length(e0);
i0 = [ibase, iimax, ibase, iimax, ibase, iimax, ibase, iimax, ibase, ibase]; lI=length(i0);

inpfreq =[0,   0,   0,    2.5,     0,    5,     0,   10,     0,     0]/4;
inpphase0 = [0, 0, 0, -pi/2, 0, -pi/2, 0, -pi/2, 0, 0];

numIPcombinations = 1; % only one row of input train
%% Computing number and sets of simulations with sufficient repeats of every combination of sigma, theta, inputs

sigmaE = sigmas; sigmaI = sigmas;
nsigs = numel(sigmaE);

ntheta = numel(thetaE);
E0 = 1:numIPcombinations; I0 = 1:numIPcombinations;
uniqE0 = E0(:); uniqI0 = I0(:);
nips = numel(uniqE0);

sigmasE = repmat(sigmaE(:)',[niters,nips*ntheta*1]);
sigmasI = repmat(sigmaI(:)',[niters,nips*ntheta*1]);
sigmaids = repmat((1:numel(sigmaE)), [niters, nips*ntheta*1]); sigmaids = sigmaids(:);
sigmas = [sigmasE(:),sigmasI(:)];
sigmasEI = [sigmaE(:),sigmaI(:)]; sigmasEI = unique(sigmasEI,'rows');
sigmaE = sigmasEI(:,1); sigmaI = sigmasEI(:,2);
selsigmaitermatrix = sigmas(:,1)==sigmaE(:)' & sigmas(:,2)==sigmaI(:)';


thetasE = repmat(thetaE(:)',[niters*nsigs,nips*1]);
thetasI = repmat(thetaI(:)',[niters*nsigs,nips*1]);
thetaids = repmat((1:numel(thetaE)), [niters*nsigs, nips*1]);
thetas = [thetasE(:),thetasI(:)];
thetasEI = [thetaE(:),thetaI(:)]; thetasEI = unique(thetasEI,'rows');
thetaE = thetasEI(:,1); thetaI = thetasEI(:,2);
selthetaitermatrix = thetas(:,1)==thetaE(:)' & thetas(:,2)==thetaI(:)';


E0 = repmat(uniqE0(:)',[niters*nsigs*ntheta,1]);
I0 = repmat(uniqI0(:)',[niters*nsigs*ntheta,1]);
stimids = repmat((1:nips),[niters*nsigs*ntheta,1]);
E0 = E0(:); I0 = I0(:);
selipitermatrix = (E0(:) == (uniqE0(:)')) & (I0(:) == (uniqI0(:)'));

nsimulations = numel(E0);


stimparamids = (stimids(:)*1 + thetaids(:)*1e3+ sigmaids(:)*1e6);

if ~exist('poplnFunc','var')
    poplnFunc = [];
end



% OU noise generation

noise_t = tspan;
noisefile = fullfile('./NoiseFiles',['passband_OUinp',num2str(theta_passband)],['newRectSinmaxip',num2str(imax),'baseip',num2str(ebase),'_', num2str(ibase),'/RectSin_',num2str(length(inpfreq)),'freqs_OU_lopas__',num2str(numel(thetas)),'theta_',num2str(niters),'.mat']);
noisefile
overWriteFlag2 = 1==0;

if ~exist(noisefile,'FILE')
    mkdir('./NoiseFiles');
    mkdir(fullfile('./NoiseFiles',['passband_OUinp',num2str(theta_passband)]));
    mkdir(fullfile('./NoiseFiles',['passband_OUinp',num2str(theta_passband)],['newRectSinmaxip',num2str(imax),'baseip',num2str(ebase),'_', num2str(ibase)]));
    overWriteFlag2 = 1==1;
    wnoise = randn([numel(E0)+numel(I0),numel(noise_t)]);

    stimid = @(t) max((1:length(diff(stimwin)))*(abs(diff(t(:)'>stimwin(:),[],1))==1), 1); %@(t) find(abs(diff(t(:)'>stimwin(:)))==1);

    wnoise = wnoise.*sigmas(:);
    meanip = [ ...
        E0(:)* ( ebase  +  (e0(stimid(noise_t))-ebase) ...
        .*abs( cos(inpphase0(stimid(noise_t))+2*pi*inpfreq(stimid(noise_t)).*(noise_t-stimwin(stimid(noise_t)))))) ;
        I0(:)* ( ibase  +  (i0(stimid(noise_t))-ibase) ...
        .*abs( cos(inpphase0(stimid(noise_t))+2*pi*inpfreq(stimid(noise_t)).*(noise_t-stimwin(stimid(noise_t)))))) ];

    wnoise = wnoise.*meanip;
    infind = isinf(thetas(:));
    thetas_woinf = thetas(:); thetas_woinf(infind) = 1;

    % setup and run OU_noise input generation model
    OUpop = OU_noise_lopass_form(nsimulations, {wnoise, meanip, noise_t}, thetas_woinf(:));
    OUpop.input([], noise_t);
    solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);

    OUpop.ode(solver);
    noise_t = OUpop.EIpairs.t;
    noise = OUpop.EIpairs.R;

    % when theta is inf, replace with white noise
    noise(infind,:) = wnoise(infind, :) + meanip(infind,:);
    clear OUpop;

    save(noisefile,'noise_t','wnoise','noise','stimparamids','niters','uniqE0','uniqI0','thetaE','thetaI','sigmaE','sigmaI','selipitermatrix','selsigmaitermatrix','selthetaitermatrix','-v7.3');
    clear noise_t wnoise noise
end

% Instantiate JS population with link to noisefile
JS_pop = ISN_KkSR_JS_OUip(nsimulations, 'normadd', noisefile);


%% Specify simulation timepoints and 0 baseline external input
JS_pop.input(0*[E0(:);I0(:)], tspan);

popln = JS_pop
clear JS_pop
%%
 %% Running simulation
    stimParamIDs = stimparamids;
    solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);

    if overWriteFlag1 || overWriteFlag2
    %         [t,y,lfp] = GammaHarmonicsTest_JSnoisy_1500ms_params_stimperiod(popln, solver, [-1 -1], stimParamIDs, niters, runflag);%save([basefolder,'/simdata/',filename,'.mat'],'-v7.3');
        if ~isempty(solver)
            if exist('odeopts','var')
                popln.ode(solver, odeopts)
            else
                popln.ode(solver);
            end
        else
            popln.ode();
        end
        % appending LFP proxy (-rE-rI) to the firing rate data
        popln.EIpairs.R = [popln.EIpairs.R; -popln.EIpairs.R(1:end/2,:)-popln.EIpairs.R(end/2+1:end,:)];
        
        mkdir([basefolder,'/simdata/']);
        save([basefolder,'/simdata/',filename,'.mat'],'-v7.3');
    else
        load([basefolder,'/simdata/',filename,'.mat']);
    end

    %% plot
    f = figure('WindowState','maximized','InvertHardcopy','on');
    subplot(3,1,1); plot(popln.EIpairs.t, popln.EIpairs.R(2*end/3+1:end,:)); xlim([0 3]); ylabel('LFP proxy (-rE-rI)');
    
    load(noisefile);
    subplot(3,1,2); plot(noise_t, noise(1:end/2,:)); xlim([0 3]); ylabel('i_E');
    subplot(3,1,3); plot(noise_t, noise(end/2+1:end,:)); xlim([0 3]); ylabel('i_I');
    xlabel('Time (s)');
    title(['OU i/p cutoffs(Hz): E: ', num2str(thetaE/2/pi), '  I: ', num2str(thetaI/2/pi)]);

    % Figure formatting for paper and save
    f= gcf; postformatFig;
    savefig([basefolder,'/lfpproxy_rectsine',num2str(thetaE/2/pi),'_',num2str(thetaI/2/pi),'.fig']);
    print(f,[basefolder,'/lfpproxy_rectsine',num2str(thetaE/2/pi),'_',num2str(thetaI/2/pi),'.png'],'-dpng','-r600');


    %% Gamma burst
    
    y = popln.EIpairs.R; t = popln.EIpairs.t;
    
    t=popln.EIpairs.t; y=popln.EIpairs.R; 
    if ~exist('doMP','var')
        doMP = 1==0;
    end

    if doMP    
        homedir = pwd();
        runtempdir = './MPruntemp_simulations_Batchiters_sigma_1ip_ou';
        mkdir(runtempdir); cd(runtempdir)
        
        % MP
        maxIteration=100;
        adaptiveDictionaryParam=0.9;
        dictionarySize=2500000;
        [gaborInfo,header] = getStochasticDictionaryMP3p1_parallel(y,t,maxIteration,adaptiveDictionaryParam,dictionarySize);
        
        mptimeVals = t;
        clear y;
        cd(homedir);
        mkdir(fullfile(basefolder,'MP'));
        gaborfile = fullfile(basefolder,'MP',[filename,'_3proxies.mat']);
        
        try
            save(gaborfile, 'gaborInfo','header','mptimeVals','niters','stimParamIDs','e0','i0','sigmas','-v7.3');
        catch
            save(gaborfile, 'gaborInfo','header','mptimeVals','niters','stimParamIDs','e0','i0','sigmas','-v7.3');
        end
                
                % load presaved gabor - multiple tries in case of exception
                % when multiple versions of the code are run in different
                % processes
                try
                    load(gaborfile, 'gaborInfo','header','niters','stimParamIDs');
                catch
                    try
                        load(gaborfile, 'gaborInfo','header','niters','stimParamIDs');
                    catch
                        try
                            load(gaborfile, 'gaborInfo','header','niters','stimParamIDs');
                        catch
                            load(gaborfile, 'gaborInfo','header','niters','stimParamIDs');
                        end
                    end
                end

    %            % gamma length analysis
    %             [lengthList,freqList,timeList,alist,trialList,~,~] = getBurstLengthMP_simulations([], t, gaborInfo, header, [10, 70]);
    
    end

end

% Function to map simulation time to input parameters
function [val] = stimid(t,win)
    val = max(1,(1:length(win)-1)*(abs(diff(t(:)'>win(:),   [], 1))==1));
end