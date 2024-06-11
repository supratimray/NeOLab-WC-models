% Simulating sustained constant input without noisy fluctuations

% poplnFunc = 'JS2014' 
%               or 'KkSR' for the JS model used in Krishnakumaran and Ray, Cerebral Cortex 2023.

% poplnFunc = model object or anonymous function instantiating a model object executable as GenerateNeuronPopulation_object = poplnFunc(nsimulation);
%%
if ~exist('modelname','var')
    modelname = [];
end
if isempty(modelname)
    modelname = 'JS2014'; % 'JS2014' or 'KkSR'
end
e0 = 0:0.5:22; lE=length(e0);
i0 = 0:0.5:22; lI=length(i0);
tspan = (0:2000)*1e-3; % seconds; Simulation timestamps

%% setting up solver
eulersolver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
% assign solver=[] for default ode45 solver in Matlab
   %  or solver=eulersolver for basic euler-forward solver 
   %  or solver=@(fn, t, y) <solver of choice>(fn, t, y,
   %  <additionaloptions>) to use custom-defined solver
solver = []; 

%% unwrapping input combinations
[E0, I0] = meshgrid(e0, i0);
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);

%% Instantiating population
if ~exist('modelname','var')
    modelname = 'JS2014';
end

if strcmp(modelname,'JS2014')
    model = @(N) ISN_JS2014(N);
elseif strcmp(modelname,'KkSR')
    model = @(N) ISN_KkSR_JS_OUip(N);
end

JS_pop = model(nsimulations);

%% savedir setup
basefolder = ['./',JS_pop.Population_Name,'_ConstantIP_results'];
filename = 'Simdata';
mkdir(basefolder);
savename = fullfile(basefolder,[filename,'.mat']);
if ~exist(savename,'file')
%% Loading constant inputs throughout tspan
JS_pop.input([E0;I0], tspan);

%% Running simulation

if exist('solver','var')
    
    if exist('odeopts','var')
       JS_pop.ode(solver, odeopts)
    else
        JS_pop.ode(solver);
    end
else
    JS_pop.ode(); % uses ode45 function
end


%% analysis
sel = JS_pop.EIpairs.t>0.25;
t = JS_pop.EIpairs.t(sel);
re = JS_pop.EIpairs.R(1:end/2,sel);
ri = JS_pop.EIpairs.R(end/2+1:end,sel);
lfp = -re -ri;

efr = mean(re,2);
ifr = mean(ri,2);

lfpfft = 2*log10(abs(fft(lfp,[],2)));
f = (0:size(lfpfft,2)-1)/(t(end)-t(1));
gammasel = f>=30 & f<=80;
[pkpower,pkfreqid] = max(lfpfft(:,gammasel),[],2);
fgammasel = f(gammasel);
pkfreq = fgammasel(pkfreqid);

%% save and plot
E0 = reshape(E0,[lI,lE]);
I0 = reshape(I0,[lI,lE]);
efr = reshape(efr,[lI,lE]);
ifr = reshape(ifr,[lI,lE]);
pkpower = reshape(pkpower,[lI,lE]);
pkfreq = reshape(pkfreq,[lI,lE]);
disp(['Saving simulation record to... ', savename]);
save(savename,'-v7.3');
else
load(savename);
end
figure('windowstate','maximized');
subplot(2,3,1); pcolor(E0,I0,efr); colormap jet; colorbar; shading interp; clim([0 1]); 
title('E firing rate'); xlabel('I_E'); ylabel('I_I');
subplot(2,3,2); pcolor(E0,I0,ifr); colormap jet; colorbar; shading interp; clim([0 1]);
title('I firing rate'); xlabel('I_E'); ylabel('I_I');

subplot(2,3,4); pcolor(E0,I0,pkpower); colormap jet; colorbar; shading interp; clim([prctile(pkpower,10,"all") prctile(pkpower,90,"all")]); 
title({'peak gamma power','(log_{10}(power))'}); xlabel('I_E'); ylabel('I_I');
subplot(2,3,5); pcolor(E0,I0,pkfreq); colormap jet; colorbar; shading interp; clim([40 70]);
title({'Peak gamma freq','(Hz)'}); xlabel('I_E'); ylabel('I_I');
%% Phase space analysis inputs 
if strcmp(modelname,'JS2014')
    IE_ToPlot = 8;
    II_ToPlot = 8;
elseif strcmp(modelname,'KkSR')
    IE_ToPlot = 8.5;
    II_ToPlot = 6.5;
end
% example trajectory to plot
selid = (E0(:)==IE_ToPlot) & (I0(:)==II_ToPlot);
rEsel = JS_pop.EIpairs.R(find(selid),:);
rIsel = JS_pop.EIpairs.R(find(selid)+end/2,:);
%% Phase diagram plotting
figure_singleInput = figure('WindowState','maximized','InvertHardcopy','on');
subplot(1,2,1);
% run describe dynamics on @(nsims) ISN_KkSR_JS_OUip(nsims)
NCvariableIDs = [1,2];
rbounds = [0, 1];
describeDynamics(figure_singleInput, gca, model, [IE_ToPlot,II_ToPlot], NCvariableIDs, {rbounds, rbounds});
hold on;
plot(rEsel(1,:),rIsel(1,:),'color',[[0.5,0.5,0.5]],'displayname','Trajectory example'); %#ok
scatter(rEsel(1,1), rIsel(1,1),[],[[0.5,0.5,0.5]],'marker','x');                        %#ok
scatter(rEsel(1,end), rIsel(1,end),[],'k','filled','marker','o');
title({[JS_pop.Population_Name,'  :  Phase diagram'], ['I_E = ',num2str(IE_ToPlot)], ['I_I = ',num2str(II_ToPlot)]}); 
xlabel('r_E'); ylabel('r_I');
legend({'Field', 'E nullcline', 'I nullcline', 'Trajectory example'});

