%%

rs = 0:1:40; rs = rs/40*5;
log10cs = (-2:2/(numel(rs)-1):0);
cs = 10.^log10cs;
e0 = cs;
i0 = rs;
lE=length(cs);
lI=length(rs);

%%

[CS0, RS0] = meshgrid(cs, rs);
crp0 = (CS0.^2)./( ((0.3).^2) + (CS0.^2) );
E0 = 40*crp0; I0 = 32*crp0;
E0 = E0(:); I0 = I0(:);
nsimulations = length(E0);

if ~exist('poplnFunc','var')
    poplnFunc = [];
end
if isempty(poplnFunc)
    JXK2013_pop = ISN_JXK(nsimulations,RS0(:));
else
    JXK2013_pop = poplnFunc(nsimulations);
end
tspan = (0:(1+ 1 )*1e4)*1e-4; % seconds
JXK2013_pop.input([E0;I0;zeros(size(E0))], tspan);
disp(JXK2013_pop)

%%
solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
[~,~,eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft, AllData] = JXK_GammaHarmonicsTest(JXK2013_pop, solver);

%%
popln = JXK2013_pop;
Rs = popln.EIpairs.R; ts = popln.EIpairs.t;
mkdir([basefolder,'/simValue_saves']);
save([basefolder,'/simValue_saves/',filename,'.mat'], 'Rs', 'ts', 'eFR','iFR','peakFreq','harmonicFreq','freqRatio','peakAmp','harmonicAmp','powerRatio','hilbPhaseDiff', 'f_fft', 'R_fft');
%%
%popln.Population_Name = [popln.Population_Name,' _ ', filename];
disp('Going into plotGammaHarmonicTestFigures.m')
plotGammaHarmonicTestFigures_cvsr;

