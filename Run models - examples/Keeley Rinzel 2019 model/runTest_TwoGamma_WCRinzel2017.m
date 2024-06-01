function [regimevals, donetime] = runTest_TwoGamma_WCRinzel2017(poplnFunc, beta, filename, basefolder)
%%
e0 = 0:0.5:20; lE=length(e0);
i0 = 0:0.5:20; lI=length(i0);

%% Model saturates above 1
e0 = e0/20;
i0 = i0/20;

%%
[E0, I0] = meshgrid(e0, i0);
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);
if ~exist('beta','var')
    beta = 0.33;
end
if ~exist('poplnFunc','var')
    poplnFunc = [];
end
if isempty(poplnFunc)
    WCRinzel17_pop = ISN_TwoGamma_WCRinzel2017(nsimulations);
else
    WCRinzel17_pop = poplnFunc(nsimulations);
end
tspan = (0: 2e4)*1e-4; % seconds

WCRinzel17_pop.input([E0;(1-beta)*I0; (beta)*I0; zeros(size(E0)); zeros(size(E0));zeros(size(I0))], tspan);
disp(WCRinzel17_pop)
%%

% solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
[~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft] = GammaHarmonicsTest(WCRinzel17_pop); %#ok

popln = WCRinzel17_pop;
popln.Population_Name = [popln.Population_Name,'_beta', num2str(beta)];
figure; plot(popln.EIpairs.t,  popln.EIpairs.R((E0(:)==5)&(I0(:)==5),:));
xlabel('time(s)'); ylabel('E pop activity');
title({'Dual gamma Keeley-Rinzel 2017 model',['Competition factor (\beta) = ', num2str(beta)]})
plotGammaHarmonicTestFigures;
%%
end