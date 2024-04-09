%% Run vanilla Jadi and Sejnowski model with 48 different constant inputs
runTest_WCJS2014; 

%% Run Krishnakumaran and Ray 2023 JS-model (slightly changed parameters) with 48 different constant inputs
% constant DC-only inputs
runTest_WCJS2014('KkSR'); 

%% Run 50 iterations of custom inputs in Krishnakumaran and Ray 2023 JS-model with noisy OU inputs
% rectified sinusoidal inputs
run_JS_OUip_rectifiedsineinput(50); 

