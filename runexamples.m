%% Run Jadi and Sejnowski model with 48 different constant inputs
modelname = 'JS2014'; 
run_JS2014_allInputCombinations; 

%% Run Krishnakumaran and Ray 2023 JS-model (slightly changed parameters) with 48 different constant inputs
% constant DC-only inputs
modelname = 'KkSR';
runTest_WCJS2014('KkSR'); 

%% Run multiple iterations of custom time-varying inputs in Krishnakumaran and Ray 2023 JS-model with noisy OU inputs
% rectified sinusoidal inputs
run_KkSR_rectsine;

%% Run phase-space analysis and parameter-sweep (for instance, toggle between JS2014 and KkSR2023 type input configurations)
run_JS_KkSR_noisy_paramsweep;