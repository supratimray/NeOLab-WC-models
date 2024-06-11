%% Run Jadi and Sejnowski model with 48 different constant inputs
modelname = 'JS2014'; 
run_JS2014_allInputCombinations; 

clear all;
%% Run Krishnakumaran and Ray 2023 JS-model (slightly changed parameters) with 48 different constant inputs
% constant DC-only inputs
modelname = 'KkSR';
run_JS2014_allInputCombinations;

clear all;
%% Run multiple iterations of custom time-varying inputs in Krishnakumaran and Ray 2023 JS-model with noisy OU inputs
% rectified sinusoidal inputs
run_KkSR_rectsine; % TODO correct noise parameters

clear all;
%% Run phase-space analysis and parameter-sweep (for instance, toggle between JS2014 and KkSR2023 type input configurations)
run_JS_KkSR_noisy_paramsweep; % TODO correct noise parameters

clear all;
%% TODO JXK, Keeley Rinzel 2017 dual gamma