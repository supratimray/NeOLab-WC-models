function runTest_WCJS2014(poplnFunc, filename, basefolder)
% Simulating sustained constant input without noisy fluctuations

% poplnFunc = 'JS2014' 
%               or 'KkSR' for the JS model used in Krishnakumaran and Ray, Cerebral Cortex 2023.

% poplnFunc = model object or anonymous function instantiating a model object executable as GenerateNeuronPopulation_object = poplnFunc(nsimulation);
%%
e0 = 0:0.5:20; lE=length(e0);
i0 = 0:0.5:20; lI=length(i0);

eFR = zeros(lE,lI);
iFR = zeros(lE,lI);
peakFreq = zeros(lE,lI);
harmonicFreq = zeros(lE,lI);
freqRatio = zeros(lE,lI);
peakAmp = zeros(lE,lI);
harmonicAmp = zeros(lE,lI);
powerRatio = zeros(lE,lI);
hilbPhaseDiff = zeros(lE,lI);

%%
[E0, I0] = meshgrid(e0, i0);
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);


if ~exist('poplnFunc','var')
    poplnFunc = [];
end
if isempty(poplnFunc)
    JS_pop = ISN_JS2014(nsimulations);
elseif isstr(poplnFunc)
    if strcmp(poplnFunc,'JS2014')
        JS_pop = ISN_JS2014(nsimulations);
    elseif strcmp(poplnFunc,'KkSR')
        JS_pop = ISN_KkSR_JS_OUip(nsimulations);
    end
else
    JS_pop = poplnFunc(nsimulations);
end

tspan = (0:2000)*1e-3; % seconds

% WC1972_pop.input([E0;I0], tspan);
% WC1972_pop
JS_pop.input([E0;I0], tspan);

%%
solver = [];      % @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);

if exist('solver','var')
    if exist('odeopts','var')
       JS_pop.ode(solver, odeopts)
    else
        JS_pop.ode(solver);
    end
else
    JS_pop.ode(); % uses ode45 function
end

if ~exist('basefolder','var')
    basefolder = ['./',JS_pop.Population_Name,'_ConstantIP_results']
end
if ~exist('filename','var')
    filename = 'Simdata'
end
mkdir(basefolder);
savename = fullfile(basefolder,[filename,'.mat']);
disp(['Saving simulation record to... ', savename]);
save(savename,'-v7.3');
end