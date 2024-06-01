% This program does the following

% 1. Computes the power spectral density of the signal 
% 2. Find the maximum power in the range specified by gammaRange. That is the peak gamma frequency
% 3. Compute the power in a band around the peak frequency specified by gammaBandWidth
% 4. Find the maximum power near 2xpeakGammaFreq. Find the peak harmonic power in the same way as step 3
% 5. Band-pass filter the signal around gamma and harmonic peak frequencies
% 6. Compute the analytic signals using hilbet tranform and get the phases



function [peakGammaFreq,harmonicFreq,freqRatio,gammaAmp,harmonicAmp,powerRatio,hilbPhaseDiff, freqVals, fftx]=getGammaAndHarmonicProperties(x,gammaRangeHz,gammaBandwidthHz,tMS,xE,xI)
    EIphaseflag = exist('xE','var') & exist('xI','var');
    if ~exist('gammaRangeHz','var');        gammaRangeHz = [30 75];         end
    if ~exist('gammaBandwidthHz','var');    gammaBandwidthHz = 20;          end

    delta = gammaBandwidthHz/2;
    BW = gammaBandwidthHz;
    
    Fs = round(1000./(tMS(2)-tMS(1)));
    T  = (tMS(end)-tMS(1))/1000; % Duration in seconds
    freqVals = (0:1:(size(x,2)-1))/T;

    fftx = fft(x,[],2);
    gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
    fftGamma = fftx(:,gammaRangePos);
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
    hilbPhaseDiff = [];
    EIphasediff = [];
    for i = 1:size(x,1)
        
        n = 4;
        [B,A] = butter(n,[peakGammaFreq(i)-delta, peakGammaFreq(i)+delta]/(Fs/2));
        [D,C] = butter(n,[harmonicFreq(i)-delta, harmonicFreq(i)+delta]/(Fs/2));
        gammaSig = filtfilt(B,A,x(i,:));
        if EIphaseflag
            gammaE = filtfilt(B,A,xE(i,:));
            gammaI = filtfilt(B,A,xI(i,:));
        end
        harmonicSig = filtfilt(D,C,x(i,:));

        G = hilbert(gammaSig);
        if EIphaseflag
            GxE = hilbert(gammaE);
            GxI = hilbert(gammaI);
            EIphasediffval= wrapTo180(rad2deg(circ_mean(angle(GxE) - angle(GxI),[],2)));
            EIphasediff = [EIphasediff;EIphasediffval];
        end
        H = hilbert(harmonicSig);
        phaseDiffFiltered = (angle(H)-2*angle(G));
        [~,rho] = rose(phaseDiffFiltered,25);
        meanPhaseDiff = getPhaseProperties(rho);
        %ph = [meanPhaseDiff-360, meanPhaseDiff];
        hilbPhaseDiff = [hilbPhaseDiff;meanPhaseDiff]; %[hilbPhaseDiff; ph(find(abs(ph)<=180,1))];
    end
    if EIphaseflag
%         EIphasediff = wrapTo180(rad2deg(circ_mean(angle(GxE) - angle(GxI),[],2)));
        fold = gcf;
        fphaseEI = figure('windowstate','maximized');
        e0 = 0:0.5:20; lE=length(e0);
        i0 = 0:0.5:20; lI=length(i0);
        [E0, I0] = meshgrid(e0, i0);
        E0 = E0(:); I0 = I0(:);
        nsimulations = length(E0);
        subplot(2,3,1); pcolor(e0,i0,reshape(EIphasediff,[lI,lE])); shading flat;
        caxis([-180,180]); colormap(fphaseEI, hsv);
        colorbar;
        filename = 'N:\Students\R. Krishnakumaran\RK NeOLab models_ correct\Testfigs';
        mkdir(filename); filenamt = [filename, '\fphaseEIfig'];
        savefig(fphaseEI, [filename,'.fig']);
        print([filename,'.tif'],'-dtiff','-r480');
        
        figure(fold);
    end
end