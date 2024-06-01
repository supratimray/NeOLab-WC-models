% modelParam: 'sig', 'pp' or 'pl' for sigmoidal, piecewise power or piecewise linear

function [x, yy, eFR,iFR,peakFreq,harmonicFreq,freqRatio,peakAmp,harmonicAmp,powerRatio,hilbPhaseDiff, f_fft, R_fft] = GammaHarmonicsTest(nrnpop,solver, odeopts,displayFlag)

if ~exist('displayFlag','var');        displayFlag=0;                   end

eqnName = nrnpop.Population_Name;
%[wcParams,stimParams] = defaultParams_JS2014; % Get default parameters
%
%stimParams.e = e0;
%stimParams.i = i0;
%wcParams.modelParam = 'sig';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodTimeVals = [1001:2000]*1e-3;    % Compute parameters for the last 1 second

%y0 = [0 0]; % start from origin
%[t,y] = ode45(@(t,y) eqn_WCJS2014(t,y,wcParams,stimParams),tVals,y0);
if exist('solver','var')
    if exist('odeopts','var')
       nrnpop.ode(solver, odeopts)
    else
        nrnpop.ode(solver);
    end
else
    solver = @(updateFn, tlist, y0) eulerMethod(updateFn, tlist, y0);
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
% if displayFlag
%     
%     figure();
%     subplot(221);
%     plot(t,y(:,1),'r'); hold on;
%     plot(t,y(:,2),'b');
%     plot(t,eFR+zeros(1,length(t)),'r--'); 
%     plot(t,iFR+zeros(1,length(t)),'b--');
%     xlabel('Time (ms)');
% 
%     fftye = fft(y(goodTimePos,(Ne-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups)));
%     fftyi = fft(y(goodTimePos,(Ni-1)*nrnpop.EIpairs.NNeuronGroups+(1:nrnpop.EIpairs.NNeuronGroups)));
%     freqVals = 0:999; % Hz
%     
%     subplot(222);
%     plot(freqVals,log10(abs(fftye)),'r'); hold on;
%     plot(freqVals,log10(abs(fftyi)),'b');
%     plot(peakFreq,log10(peakAmp),'ro');
%     plot(harmonicFreq,log10(harmonicAmp),'bo');
%     xlim([0 150]);
%     xlabel('Frequency (Hz)');
%     ylabel('log Power');
%     
%     subplot(212);
%     xRange = [0 1];
%     yRange = [0 1];
%     numDivisions = 20;
%     xList = linspace(xRange(1),xRange(2),numDivisions+1);
%     yList = linspace(yRange(1),yRange(2),numDivisions+1);
%     
%     describeDynamics(eqnName,xList,yList,wcParams,stimParams);
% end
end