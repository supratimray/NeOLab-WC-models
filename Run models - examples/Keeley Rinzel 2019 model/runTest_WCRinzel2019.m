clear;

e0 = 0:0.5:20; lE=length(e0);
i0 = 0:0.5:20; lI=length(i0);


%%
[E0, I0] = meshgrid(e0, i0);
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);
% WC1972_pop = ISN_WC1972(nsimulations);

if ~exist('poplnFunc','var')
    WCRinzel19_pop = ISN_WCRinzel2019(nsimulations);
else
    WCRinzel19_pop = poplnFunc(nsimulations);
end
tspan = [0, 2000]*1e-3; % seconds

WCRinzel19_pop.input([zeros(size(E0));zeros(size(I0));E0;I0], tspan);
disp(WCRinzel19_pop)
%%
[~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft] = GammaHarmonicsTest(WCRinzel19_pop);

goodHarmonicPos = find(harmonicAmp > 1.5);
goodGammaPos = find(peakAmp > 19);
goodPos = intersect(goodHarmonicPos,goodGammaPos);

popln = WCRinzel19_pop;
if ~exist('filename','var')
    filename = '';
end
popln.Population_Name = [popln.Population_Name,' _ ', filename];
plotGammaHarmonicTestFigures;
