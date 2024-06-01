% modelParam: 'sig', 'pp' or 'pl' for sigmoidal, piecewise power or piecewise linear

function [x, yy, eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft] = GammaHarmonicsTest(nrnpop,solver, odeopts)

if ~exist('displayFlag','var');        displayFlag=0;                   end

eqnName = nrnpop.Population_Name;
%[wcParams,stimParams] = defaultParams_JS2014; % Get default parameters
%
%stimParams.e = e0;
%stimParams.i = i0;
%wcParams.modelParam = 'sig';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodTimeVals = (1001:2000)*1e-3;    % Compute parameters for the last 1 second

%y0 = [0 0]; % start from origin
%[t,y] = ode45(@(t,y) eqn_WCJS2014(t,y,wcParams,stimParams),tVals,y0);
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
% x = y((Ne-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups), goodTimePos);
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
xE = y((Ne-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups), goodTimePos); 
xI = y((Ni-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups), goodTimePos);
xE = xE(:,samples); xI = xI(:,samples);
tMS = tgood(samples)*1e3; % converting seconds to milliseconds

gammaRangeHz = [30 75];
%%
[peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff,f_fft,R_fft]=getGammaAndHarmonicProperties(x,gammaRangeHz,20,tMS,xE,xI);

end