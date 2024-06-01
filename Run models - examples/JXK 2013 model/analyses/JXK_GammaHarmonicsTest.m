% modelParam: 'sig', 'pp' or 'pl' for sigmoidal, piecewise power or piecewise linear

function [x, yy, eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft, AllData] = JXK_GammaHarmonicsTest(nrnpop,solver, odeopts,displayFlag)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodTimeVals = (10001:20000)*1e-4;    % Compute parameters for the last 1 second
if exist('solver','var')
    if exist('odeopts','var')
       nrnpop.ode(solver, odeopts)
    else
        nrnpop.ode(solver);
    end
else
%     solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
    solver = @(updateFn, tlist, y0) ode113(updateFn, tlist, y0);
    nrnpop.ode(solver);
end

t = nrnpop.EIpairs.t; 
y = nrnpop.EIpairs.R;
goodTimePos = (t>=goodTimeVals(1)) & (t<=goodTimeVals(end));

% Get FFT of last 1 second. We get 1 second resolution
Ne = 1; Ni = 2; % The state variable in a NeuronGroup corresponding to Ecell and Icell firing rates. Ne = 1; Ni = 2; for EI pair out of NStateVarsPerGroup=2
eFR = mean(y(...
    (Ne-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups),...
    goodTimePos),...
    2);
iFR = mean(y(...
    (Ni-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups),...
    goodTimePos),...
    2);

yy = y(:, goodTimePos);
%% computing LFP proxy
x = y((Ne-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups), goodTimePos) ...
    + y((Ni-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups),goodTimePos);
x = -x;
%%
tgood = t(goodTimePos);
samples = [];
for val = goodTimeVals
    samples = [samples,find( (abs(tgood-val) == min(abs(tgood-val))) ,1)];
end
x = x(:,samples);
tMS = tgood(samples)*1e3; % converting seconds to milliseconds

gammaRangeHz = [30 75];
%%
[peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff,f_fft,R_fft, AllData]=JXK_getGammaAndHarmonicProperties(x,gammaRangeHz,20,tMS);

end