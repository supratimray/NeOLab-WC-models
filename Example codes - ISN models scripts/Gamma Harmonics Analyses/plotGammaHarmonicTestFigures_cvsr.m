if ~exist('basefolder','var')
    basefolder = './plotGammaHarmonicTestFigures_DefaultBasefolder';
end
donetime = datevec(now());
if ~exist('filename','var')
    filename = [];
end
if isempty(filename)
    filename = datestr(donetime, 'dd-mm-yyyy HH-MM-SS');
end

%% time unit correction
tstep = diff(popln.EIpairs.t(1:2));
peakAmp = peakAmp*tstep;
harmonicAmp = harmonicAmp*tstep;
R_fft = R_fft*tstep;

%%
basefilename = filename;
Resultsdir = basefolder;%[basefolder,'/Simulation Outputs'];
mkdir(Resultsdir);
savedir = [Resultsdir, '/Plots'];
mkdir(savedir);
plotsavedir = savedir;%[savedir,'/',popln.Population_Name];
mkdir(plotsavedir);
plotsavename = basefilename;

savedir = [Resultsdir, '/Simulated Data'];
mkdir(savedir);
Objsavedir = savedir;%[savedir,'/',popln.Population_Name];
mkdir(Objsavedir);
Objsavename = basefilename;

if exist('AllData','var')
    disp('AllData found')
    save([Objsavedir,'/','AllData_',Objsavename,'.mat'], '-struct', 'AllData');
end
%%
if ~exist('dontPlot','var')
    dontPlot = false;
elseif ~isbool(dontPlot)
    dontPlot = false;
end
%%
if ~dontPlot
    %% selecting regime
    % goodHarmonicPos = find(harmonicAmp > 1.5);
    % goodGammaPos = find(peakAmp > 19);
    % goodPos = intersect(goodHarmonicPos,goodGammaPos);
    goodHarmonicPos = find(harmonicAmp > 1e-6);
    goodGammaPos = find(peakAmp > 1e-3);
    goodPos = intersect(goodHarmonicPos,goodGammaPos);
    
    
    freqRatioSelected = zeros(lI,lE)/0;
    powerRatioSelected = zeros(lI,lE)/0;
    hilbPhaseDiffSelected = zeros(lI,lE)/0;
    freqRatioSelected(goodPos) = freqRatio(goodPos);
    powerRatioSelected(goodPos) = powerRatio(goodPos);
    index = hilbPhaseDiff<0; hilbPhaseDiff(index) = hilbPhaseDiff(index)+360;
    hilbPhaseDiffSelected(goodPos) = hilbPhaseDiff(goodPos);
    %
    % goodGHPos = zeros(length(goodPos),2);
    % for i = 1:length(goodPos)
    %     [row, col] = find(harmonicAmp == harmonicAmp(goodPos(i)));
    %     goodGHPos(i,1) = row;
    %     goodGHPos(i,2) = col;
    % end
    
    ingoodpos = zeros([lI,lE]);
    ingoodpos(goodPos) = 1;
    interval = 45/2; %degrees
    zeroPhaseDiffregion = (reshape(hilbPhaseDiffSelected,[lI, lE]) <= 180+ interval) & (reshape(hilbPhaseDiffSelected,[lI, lE]) >= 180-interval);
    
    %% Firing Rate plot
    
    figure('name',popln.Population_Name);
    set(gcf, 'Windowstate', 'normal');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    colormap jet
    E0 = (reshape(E0, [lI, lE])); I0 = (reshape(I0, [lI, lE]));
    subplot(121)
    pcolor(e0,i0,(reshape(eFR,[lI, lE]))); shading interp;
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('E firing rate');
    %xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    xlabel('Stimulus contrast (c)','FontWeight','bold'); ylabel('Stimulus size (r)','FontWeight','bold');

    subplot(122)
    pcolor(e0,i0,(reshape(iFR,[lI, lE]))); shading interp;
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('I Firing rate');
    xlabel('Stimulus contrast (c)','FontWeight','bold'); %ylabel('Stimulus size (r)','FontWeight','bold');
    % xlabel('Input to E pop. (i_{E})','FontWeight','bold');
    %%
    print([plotsavedir, '/', 'FiringRates_', plotsavename,'.png'],'-dpng','-r480');
    savefig(gcf, [plotsavedir, '/', 'FiringRates_', plotsavename,'.fig'] ,'compact');
    % print([plotsavedir, '/', 'FiringRates_', plotsavename,'.svg'],'-dsvg','-r480');
    % print([plotsavedir, '/', 'FiringRates_', plotsavename,'.tif'],'-dtiff','-r480');
    
    %% Gamma and Harmonic Analysis plot
    figure('name',popln.Population_Name);
    set(gcf, 'Windowstate', 'maximize');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    colormap jet
    subplot(231)
    pcolor(e0,i0,(reshape(peakAmp,[lI, lE]))); shading interp;
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('Peak Gamma Amplitude');
    
    subplot(234)
    pcolor(e0,i0,(reshape(harmonicAmp,[lI, lE]))); shading interp;
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    xlabel('Stimulus contrast (c)','FontWeight','bold'); ylabel('Stimulus size (r)','FontWeight','bold');
    % xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    title('Harmonic Amplitude');
    
    subplot(232)
    pcolor(e0,i0,(reshape(peakFreq,[lI, lE]))); shading interp;
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('Peak Frequency (Hz)');
    
    subplot(235)
    pcolor(e0,i0,(reshape(freqRatioSelected,[lI, lE]))); shading interp;
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    caxis([0 3])
    % hold on;
    % [~,contourline] = contour(e0,i0, ~isnan(reshape(freqRatioSelected,[lI, lE])), 1, 'c');
    % contourline.LineWidth = 2;
    title('Frequency Ratio (H/F)');
    
    subplot(233)
    pcolor(e0,i0,(reshape(powerRatioSelected,[lI, lE]))); shading interp;
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    % xlabel('Stimulus contrast (c)','FontWeight','bold'); ylabel('Stimulus size (r)','FontWeight','bold');
    colorbar;
    % hold on;
    % [~,contourline] = contour(e0,i0, ingoodpos, 1, 'c');
    % contourline.LineWidth = 2;
    title('Power Ratio (H/F)');
    %%
    subplot(236)
    colormap(gca, hsv);
    pcolorSurface = pcolor(e0,i0,(reshape(wrapTo180(hilbPhaseDiffSelected),[lI, lE]))); % pcolorSurface.FaceColor = 'interp';
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    colormap(gca, hsv);
    colorbar;
    caxis([-180 180]);
    %%
    res = 5; %bits/color
    x = 0:2^(3*res)-1; x=x(:);
    rgbmap = [floor(x/2^(2*res))/(2^res-1), floor((x-floor(x/2^(2*res))*2^(2*res))/2^(1*res))/(2^res-1), floor((x-floor(x/2^(1*res))*2^(1*res))/2^0)/(2^res-1)];
    hsvmap = hsv(numel(x));
    noNanIndices = ~isnan(hilbPhaseDiffSelected)==1;
    hsvdata = zeros(numel(hilbPhaseDiffSelected),3)/0;
    hsvdata(noNanIndices,:) = hsvmap(round((wrapTo180(hilbPhaseDiffSelected(noNanIndices))+180)/360*(size(hsvmap,1)-1)+1),:);
    hexdata = ( hsvdata * (2.^(res*[2;1;0])) ) * (2^res-1);
    % colormap(gca, rgbmap);
    pcolorSurface = pcolor(e0, i0, reshape( hexdata(:), [lI,lE])); shading interp;
    pcolorSurface.FaceColor = 'interp';
    colormap(gca, rgbmap);
%     hold on;
    % [~,contourline] = contour(e0,i0, ingoodpos & ~(zeroPhaseDiffregion), 1, 'r');
    % contourline.LineWidth = 3;
    hold off;
    subplot(236)
    hsvdata = hsvmap(round((wrapTo180(hilbPhaseDiff)+180)/360*(size(hsvmap,1)-1)+1),:);
    hsvimg = reshape(hsvdata, [lI lE 3]);
    scale = 5;
    scale = 10; I = image(imresize(hsvimg, scale, 'bilinear'));
    linearindices = transpose(round(1:(lI-1)/(lI*scale-1):lI))*ones(1,lE*scale) + ones(lI*scale,1)*round(0:(lE-1)/(lE*scale-1):lE-1)*lI;
    noNanIndices = ~isnan(hilbPhaseDiffSelected(linearindices));
    ax = gca; colordata = ax.Children.CData;
    mesh(...
        e0(1):(e0(end)-e0(1))/(size(noNanIndices,2)-1):e0(end), ...
        i0(1):(i0(end)-i0(1))/(size(noNanIndices,1)-1):i0(end), ...
        zeros(size(noNanIndices)),...
        colordata.*noNanIndices + 0.5*(~noNanIndices)); view([0 90]);
    colormap(gca, hsv); cbar = colorbar;
%     for i = 1:numel(cbar.TickLabels)
%         cbar.TickLabels(i) = {num2str(360/(numel(cbar.TickLabels)-1)*(i-1))};
%     end
    c = colorbar;
    c.Ticks = (0:6)/6; % -180:60:180; 
    c.TickLabels = {'-180','-120',' -60','  0','  60',' 120',' 180'};
    hold on;    
    [~,contourline] = contour(e0,i0, ingoodpos & (zeroPhaseDiffregion), 1, 'k');
    contourline.LineWidth = 1.5;
    hold on;
    title('Phase Diff (H/F)');
    %%
    print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.png'],'-dpng','-r480');
    savefig(gcf, [plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.fig'] ,'compact');
    % print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.svg'],'-dsvg','-r480');
    % print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.tif'],'-dtiff','-r480');
    disp(['Plots Saved to: ', what(plotsavedir).path,'\*_',plotsavename]);
    
    %%
    InpE_regime = E0(ingoodpos & (zeroPhaseDiffregion)); InpE_regime = InpE_regime(:);
    InpI_regime = I0(ingoodpos & (zeroPhaseDiffregion)); InpI_regime = InpI_regime(:);
    Phasediff_regime = wrapTo180(hilbPhaseDiffSelected(ingoodpos & (zeroPhaseDiffregion)));
    powerRatio_regime = reshape(powerRatioSelected,[lI, lE]); powerRatio_regime = powerRatio_regime(ingoodpos & (zeroPhaseDiffregion));
    regimevals = table(sortrows([InpE_regime, InpI_regime, Phasediff_regime, powerRatio_regime], 3)) %, 'VariableNames', {'Input to E','Input to I','PhaseDifference (H-G)', 'PowerRatio (H/G)'})
    save([Objsavedir,'/','Regimevals_',Objsavename,'.mat'], 'regimevals');
    %% plot firing rates for input cases that fall within the regime
    if isempty(InpE_regime)
    else
        randind = randperm(length(InpE_regime), min(3,length(InpE_regime)));
        E_selectedcase = InpE_regime(randind); % Change this and execute section for visualizing other cases manually
        I_selectedcase = InpI_regime(randind);
        table([E_selectedcase(:), I_selectedcase(:)])%, 'VariableNames', {'Input to E','Input to I'})
        t = popln.EIpairs.t;
        timerange = t>=1.1 & t<=1.2;
        selectedcases = [];
        for i = 1:length(E_selectedcase)
            selectedcases = [selectedcases,find(E0(:)==E_selectedcase(i) & I0(:)==I_selectedcase(i))];
        end
        
        f = figure('name', 'Example Firing rates inside the estimated regime');
        set(f, 'DefaultLineLinewidth', 2.5);
        t = t(timerange);
        subplot(4,1,1); plot(t, cos(2*pi*30*t)+ 1/5*cos(2*pi*60*t)); title('Desired waveshape: cos(2{\pi} 30t) + 1/5 cos(2{\pi} 60t)');
        eR = popln.EIpairs.R(selectedcases,timerange);
        subplot(4,1,2); plot(t, eR); title('E cell firing rates')
        iR = popln.EIpairs.R(popln.EIpairs.NNeuronGroups+(selectedcases),timerange);
        subplot(4,1,3); plot(t, iR); title('I cell firing rates')
        subplot(4,1,4); plot(t, -(eR+iR)); title('E+I cell firing rates')
        xlabel('time (seconds)')
        
        %% save above example plots
        print([plotsavedir, '/', 'ExampleFiringRates_', plotsavename,'.png'],'-dpng','-r480');
        savefig(gcf,[plotsavedir, '/', 'ExampleFiringRates_', plotsavename, '.fig'],'compact');
        % print([plotsavedir, '/', 'ExampleFiringRates_', plotsavename,'.svg'],'-dsvg','-r480');
        % print([plotsavedir, '/', 'ExampleFiringRates_', plotsavename,'.tif'],'-dtiff','-r480');
        %%
    end
    
%     %% plot E firing rates in a grid for visual inspection
%     if ~exist('nr','var')
%         nr = 4;
%     end
%     if ~exist('nc','var')
%         nc = 4;
%     end
%     figure('name','E firing rates');
%     set(gcf, 'Windowstate', 'maximize');
%     set(gcf, 'color', 'w');
%     set(gcf,'InvertHardcopy','off');
%     subplot(1,1,1);
%     set(gca,'color','w');
%     plot(0,0,'w');
%     xlabel(['Subplots arranged by Stimulus contrast. (c): left-to-right =',num2str(e0(1)),' to ', num2str(e0(end))],'FontWeight','bold');
%     ylabel(['Subplots arranged by Stimulus size. (r): left-to-right =',num2str(i0(1)),' to ', num2str(i0(end))],'FontWeight','bold');
%     title('E+I cell firing rates');
%     set(gcf, 'Defaultlinelinewidth', 2);
%     
%     E0ref = E0(:); I0ref = I0(:);
%     plotrange = (popln.EIpairs.t>1 & popln.EIpairs.t <= 1+3*1/25); % t0 < t <= t0+n*1/f ==> n full cycles of an f Hz oscillation
%     subsamplestep = 4;
%     plotaxes = [];
%     
%     Ne = 1; Ni = 2; Nnrngrps = popln.EIpairs.NNeuronGroups;
%     REplusI = popln.EIpairs.R((Ne-1)*Nnrngrps + (1:Nnrngrps), plotrange) + popln.EIpairs.R((Ni-1)*Nnrngrps + (1:Nnrngrps), plotrange);
%     
%     for inde = 1:subsamplestep:length(e0)% 0:nc/length(e0):nc-1/length(e0)
%         for indi = 1:subsamplestep:length(i0)%0:nr/length(i0):nr-1/length(i0)
%             indr = indi/length(i0)*(nr-1/length(i0)); indc = inde/length(e0)*(nc-1/length(e0));
%             indr = ceil(indr); indc = ceil(indc);
%             %         [indi, inde, indr, indc, (nr-indr)*nc + indc]
%             ax = subplot(nr,nc, (nr-indr)*nc + indc);
%             plotaxes = [plotaxes, ax];
%             set(gca, 'color', [0.1, 0.1, 0.1]);
%             hold on;
%             plot(popln.EIpairs.t(plotrange), REplusI(E0==e0(inde)&I0==i0(indi),:));
%             axis tight;
%             aylim = ylim;
%             ylim([-0.05, max(aylim(2), 0.5)]);
%             hold off;
%             xlabel('time (s)');
%         end
%     end
%     linkaxes(plotaxes,'x');
%     
%     print([plotsavedir, '/', 'FiringRate grid_', plotsavename,'.png'],'-dpng','-r480');
%     savefig(gcf,[plotsavedir, '/', 'FiringRate grid_', plotsavename, '.fig'],'compact');
%     %% plot fft of E firing rates in a grid for visual inspection
%     if ~exist('nr','var')
%         nr = 5;
%     end
%     if ~exist('nc','var')
%         nc = 5;
%     end
%     figure('name','FFT of E+I firing rates');
%     set(gcf, 'Windowstate', 'maximize');
%     set(gcf, 'color', 'w');
%     set(gcf,'InvertHardcopy','off');
%     subplot(1,1,1);
%     set(gca,'color','w');
%     plot(0,0,'w');
%     xlabel(['Subplots arranged by Stimulus contrast. (c): left-to-right =',num2str(e0(1)),' to ', num2str(e0(end))],'FontWeight','bold');
%     ylabel(['Subplots arranged by Stimulus size. (r): left-to-right =',num2str(i0(1)),' to ', num2str(i0(end))],'FontWeight','bold');
%     title('FFT of E+I cell firing rates');
%     set(gcf, 'Defaultlinelinewidth', 1.5);
%     
%     E0ref = E0(:); I0ref = I0(:);
%     plotrange = (f_fft>25 & f_fft <= 150); % t0 < t <= t0+n*1/f ==> n full cycles of an f Hz oscillation
%     subsamplestep = 4;
%     plotaxes = [];
%     for inde = 1:subsamplestep:length(e0)% 0:nc/length(e0):nc-1/length(e0)
%         for indi = 1:subsamplestep:length(i0)%0:nr/length(i0):nr-1/length(i0)
%             indr = indi/length(i0)*(nr-1/length(i0)); indc = inde/length(e0)*(nc-1/length(e0));
%             indr = ceil(indr); indc = ceil(indc);
%             %         [indi, inde, indr, indc, (nr-indr)*nc + indc]
%             ax = subplot(nr,nc, (nr-indr)*nc + indc);
%             plotaxes = [plotaxes, ax];
%             set(gca, 'color', [0.1, 0.1, 0.1]);
%             hold on;
%             plot(f_fft(plotrange), abs(R_fft(E0==e0(inde)&I0==i0(indi), plotrange)));
%             axis tight;
%             aylim = ylim;
%             ylim([-0.5, max(aylim(2), 1)]);
%             hold off;
%             xlabel('frequency (Hz)');
%         end
%     end
%     linkaxes(plotaxes, 'x');
%     
%     print([plotsavedir, '/', 'FFT grid_', plotsavename,'.png'],'-dpng','-r480');
%     savefig(gcf,[plotsavedir, '/', 'FFT grid_', plotsavename, '.fig'],'compact');
    %%
    Population_Name = popln.Population_Name;
    R = popln.EIpairs.R;
    t = popln.EIpairs.t;
    W = popln.EIpairs.W;
    tau = popln.EIpairs.tau;
    Input = popln.EIpairs.Input;
    DynamicInputParameters = popln.DynamicInputParameters;
    save([Objsavedir,'/','PoplnObject_',Objsavename,'.mat'], 'Population_Name','R', 't', 'W', 'tau', 'Input', 'DynamicInputParameters','-v7.3');
end