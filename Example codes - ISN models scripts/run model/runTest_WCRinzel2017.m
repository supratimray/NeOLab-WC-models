function [regimevals, donetime] = runTest_WCRinzel2017(poplnFunc, filename, basefolder)
%%
e0 = 0:0.5:20; lE=length(e0);
i0 = 0:0.5:20; lI=length(i0);

%% Model saturates above 1
e0 = e0/20;
i0 = i0/20;
%%
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
% WC1972_pop = ISN_WC1972(nsimulations);

if ~exist('poplnFunc','var')
    poplnFunc = [];
end
if isempty(poplnFunc)
    WCRinzel17_pop = ISN_WCRinzel2017(nsimulations);
else
    WCRinzel17_pop = poplnFunc(nsimulations);
end
tspan = (0: 2e4)*1e-4; % seconds
% WC1972_pop.input([E0;I0], tspan);
% WC1972_pop
WCRinzel17_pop.input([E0;I0; zeros(size(E0));zeros(size(I0))], tspan);
WCRinzel17_pop
%%
% [~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff] = GammaHarmonicsTest(WC1972_pop, 0);
% [~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff] = GammaHarmonicsTest(pop, 0);
[~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft] = GammaHarmonicsTest(WCRinzel17_pop);

% goodHarmonicPos = find(harmonicAmp > 1.5);
% goodGammaPos = find(peakAmp > 19);
% goodPos = intersect(goodHarmonicPos,goodGammaPos);

popln = WCRinzel17_pop;
%popln.Population_Name = [popln.Population_Name,' _ ', filename];
plotGammaHarmonicTestFigures;
%%
end