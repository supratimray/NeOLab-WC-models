% This program does the following

% 1. Computes the power spectral density of the signal 
% 2. Find the maximum power in the range specified by gammaRange. That is the peak gamma frequency
% 3. Compute the power in a band around the peak frequency specified by gammaBandWidth
% 4. Find the maximum power near 2xpeakGammaFreq. Find the peak harmonic power in the same way as step 3
% 5. Band-pass filter the signal around gamma and harmonic peak frequencies
% 6. Compute the analytic signals using hilbet tranform and get the phases


function [peakGammaFreq,harmonicFreq,freqRatio,gammaAmp,harmonicAmp,powerRatio,hilbPhaseDiff, freqVals, avgfftx, AllData]=JXK_getGammaAndHarmonicProperties_paperversion(x,avgfftx,gammaRangeHz,gammaBandwidthHz,tMS)
    
    if ~exist('gammaRangeHz','var');        gammaRangeHz = [30 75];         end
    if ~exist('gammaBandwidthHz','var');    gammaBandwidthHz = 20;          end

    delta = gammaBandwidthHz/2;
    BW = gammaBandwidthHz;

    Fs = round(1000./(tMS(2)-tMS(1)));
    T  = (tMS(end)-tMS(1))/1000; % Duration in seconds
    freqVals = (0:1:(size(x,2)-1))/T;
    if isempty(avgfftx)
        fftx = fft(x,[],2);
        analysisflag = false; % no more analysis to be performed
        peakGammaFreq = [];
        harmonicFreq = [];
        freqRatio = [];
        gammaAmp = [];
        harmonicAmp = [];
        powerRatio = [];
        hilbPhaseDiff = []; 
    else
        % do gamma-harmonic extraction and phase analysis using fftx in
        % arguments
        fftx = avgfftx;
        analysisflag = true;
    end
    gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
    fftGamma = fftx(:,gammaRangePos);
    
    if analysisflag
        gammaAmp = max(abs(fftGamma), [], 2);
        gammaPos = [];
        peakGammaFreq = [];
        for i =1:size(x,1)
            gammaPos = [gammaPos;gammaRangePos(find(abs(fftGamma(i,:))==gammaAmp(i),1))];
            peakGammaFreq = [peakGammaFreq; freqVals(gammaPos(end))];
        end
        %     if isrow(peakGammaFreq)
        %         peakGammaFreq = transpose(peakGammaFreq);
        %     end
        %     gammaPhase = angle(fftx(gammaPos));
        
        estHarmonicFreq = 2*peakGammaFreq;
        harmonicPos = []; harmonicAmp = [];
        for i =1:size(x,1)
            harmonicRangePos = intersect(find(freqVals>=estHarmonicFreq(i)-BW),find(freqVals<estHarmonicFreq(i)+BW));
            fftHarmonic = fftx(i, harmonicRangePos);
            harmonicAmp = [harmonicAmp;max(abs(fftHarmonic), [], 2)];
            harmonicPos = [harmonicPos;harmonicRangePos(find(abs(fftHarmonic)==harmonicAmp(i),1))];
        end
        %    harmonicPos = harmonicRangePos(find(abs(fftHarmonic)==harmonicAmp,1));
        harmonicFreq = freqVals(harmonicPos);
        if isrow(harmonicFreq)
            harmonicFreq = transpose(harmonicFreq);
        end
        %     harmonicPhase = angle(fftx(harmonicPos));
        
        freqRatio = harmonicFreq./peakGammaFreq;
        powerRatio = (harmonicAmp./gammaAmp).^2;
        %     phaseDiff  = harmonicPhase-2*gammaPhase;
        %     hilbPhaseDiff = 180*phaseDiff;
        
        %     hilbx = hilbert(x);
        %     gammaPhaseHilb = unwrap(angle(hilbx(gammaPos)));
        %     harmonicPhaseHilb = unwrap(angle(hilbx(harmonicPos)));
        %     hilbPhaseDiff = 180*(gammaPhaseHilb-harmonicPhaseHilb);
        hilbPhaseDiff = []; PhaseDifferences_t = []; GammaAmp_t = [];
        for i = 1:size(x,1)
            
            n = 4;
            [B,A] = butter(n,[peakGammaFreq(i)-delta, peakGammaFreq(i)+delta]/(Fs/2));
            [D,C] = butter(n,[harmonicFreq(i)-delta, harmonicFreq(i)+delta]/(Fs/2));
            gammaSig = filtfilt(B,A,x(i,:));
            harmonicSig = filtfilt(D,C,x(i,:));
            
            G = hilbert(gammaSig);
            H = hilbert(harmonicSig);
            phaseDiffFiltered = (angle(H)-2*angle(G));
            gammaamplitude = abs(G).^2;
            [~,rho] = rose(phaseDiffFiltered,25);
            meanPhaseDiff = getPhaseProperties(rho);
            %ph = [meanPhaseDiff-360, meanPhaseDiff];
            hilbPhaseDiff = [hilbPhaseDiff;meanPhaseDiff]; %[hilbPhaseDiff; ph(find(abs(ph)<=180,1))];
            PhaseDifferences_t = [PhaseDifferences_t, phaseDiffFiltered(:)];
            GammaAmp_t = [GammaAmp_t, gammaamplitude(:)];
        end
    end
    
    tstep = (tMS(2)-tMS(1))/1000; % least time-step in seconds;
    if analysisflag
        AllData.peakGammaFreq = peakGammaFreq;
        AllData.harmonicFreq = harmonicFreq;
        AllData.freqRatio = freqRatio;
        AllData.gammaAmp = gammaAmp*tstep;
        AllData.harmonicAmp = harmonicAmp*tstep;
        AllData.powerRatio = powerRatio;
        AllData.hilbPhaseDiff = hilbPhaseDiff;
        AllData.hilbPhaseDiff(hilbPhaseDiff<0) = hilbPhaseDiff(hilbPhaseDiff<0)+360;
        AllData.PhaseDiff_t = PhaseDifferences_t;
        AllData.GammaAmp_t = GammaAmp_t;
    end
    AllData.freqVals = freqVals;
    AllData.fftx = fftx*tstep;
end

