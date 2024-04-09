function [regimevals, donetime] = runTest_JXK2013_JXKregime_rvsc(poplnFunc, filename, basefolder)
%%

rs = 0:1:40; rs = rs/40*5;
log10cs = (-2:2/(numel(rs)-1):0);
cs = 10.^log10cs;
e0 = cs; % e0_actual = 40*(cs.^2)/( ((0.3).^2) + (cs.^2) );
i0 = rs; % i0_actual = 32*(cs.^2)/( ((0.3).^2) + (cs.^2) );
lE=length(cs);
lI=length(rs);

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
% [E0, I0] = meshgrid(e0, i0);
[CS0, RS0] = meshgrid(cs, rs); 
crp0 = (CS0.^2)./( ((0.3).^2) + (CS0.^2) ); E0 = 40*crp0; I0 = 32*crp0;
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);

if ~exist('poplnFunc','var')
    poplnFunc = [];
end
if isempty(poplnFunc)
    JXK2013_pop = ISN_JXK(nsimulations);
else
    JXK2013_pop = poplnFunc(nsimulations);
end
tspan = (0:(1+ 1 )*1e4)*1e-4; % seconds
% tspan = (0:(1+ 10 )*1e4)*1e-4; % seconds
% tspan = (0:(1+ 30 )*1e4)*1e-4; % seconds
JXK2013_pop.input([E0;I0;zeros(size(E0))], tspan);
JXK2013_pop
% JS_pop.input([E0;I0], tspan);
% JS_pop
%%
% [~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff] = GammaHarmonicsTest(pop, 0);
% odeopts = odeset('RelTol',1e-3,'AbsTol',1e-1);
solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
% [~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft, AllData] = Modifying_GammaHarmonicsTest(JXK2013_pop, solver);
[~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft, AllData] = JXK_GammaHarmonicsTest(JXK2013_pop, solver);
% [~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff] = GammaHarmonicsTest(JS_pop, 0);
%%
popln = JXK2013_pop;
Rs = popln.EIpairs.R; ts = popln.EIpairs.t;
mkdir([basefolder,'/simValue_saves']);
save([basefolder,'/simValue_saves/',filename,'.mat'], 'Rs', 'ts', 'eFR','iFR','peakFreq','harmonicFreq','freqRatio','peakAmp','harmonicAmp','powerRatio','hilbPhaseDiff', 'f_fft', 'R_fft');
%%
% load 'tempFiring rate values etc.mat'
%popln.Population_Name = [popln.Population_Name,' _ ', filename];
disp('Going into plotGammaHarmonicTestFigures.m')
plotGammaHarmonicTestFigures_cvsr;

end