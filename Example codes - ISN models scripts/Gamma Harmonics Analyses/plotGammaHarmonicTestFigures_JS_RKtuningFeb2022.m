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
%%
Rs = popln.EIpairs.R; ts = popln.EIpairs.t;
mkdir([basefolder,'/simValue_saves']);
save([basefolder,'/simValue_saves/',filename,'.mat'], 'Rs', 'ts', 'eFR','iFR','peakFreq','harmonicFreq','freqRatio','peakAmp','harmonicAmp','powerRatio','hilbPhaseDiff', 'f_fft', 'R_fft');


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
Se = reshape(eFR.*(1-eFR), [lI, lE]); Si = reshape(iFR.*(1-iFR), [lI, lE]);
%% JS regimes
dSe_x_e0res_t = diff(Se,1,2)./diff(reshape(eFR,[lI,lE]),1,2); dSe_x_e0res = dSe_x_e0res_t(1:end-1,:);
dSi_x_i0res_t = diff(Si,1,1)./diff(reshape(iFR,[lI,lE]),1,1); dSi_x_i0res = dSi_x_i0res_t(:,1:end-1);
isuperlinear = dSi_x_i0res >0; 
esublinear = dSe_x_e0res <=0; % eweaksuperlinear = dSe_x_e0res <dSi_x_i0res;
e0_S = e0(1:end-1) + (e0(2)-e0(1))/2;
i0_S = i0(1:end-1) + (i0(2)-i0(1))/2;
dim = popln.NStateVars;
% Wee = popln.EIpairs.W(1:end/dim            ,1:end/dim);     Wei = popln.EIpairs.W(1:end/dim            ,end/dim+1:end*(2/dim));
% Wie = popln.EIpairs.W(end/dim+1:end*(2/dim),1:end/dim);     Wii = popln.EIpairs.W(end/dim+1:end*(2/dim),end/dim+1:end*(2/dim));
Wee = popln.connectivityMaxW(1,1);     Wei = popln.connectivityMaxW(1,2);
Wie = popln.connectivityMaxW(2,1);     Wii = popln.connectivityMaxW(2,2);
Trace = ( (-1+Wee*Se(:))/popln.EIpairs.tau(1) + (-1+Wii*Si(:))./popln.EIpairs.tau(end/dim+1) );
% AlternativeJSRegime = diff(reshape(Trace,[lI,lE]),1,1)>0;
TraceCondition = reshape(Trace>0, [lI,lE]);
det = (-1+Wee*Se(:)) .* (-1+Wii*Si(:)) - (Wei*Se(:)) .* (Wie*Si(:));
detCondition = reshape(det<=(Trace/2).^2, [lI,lE]); %reshape(det>0, [lI,lE]);
supHopf = TraceCondition & detCondition;
% JSRegime = AlternativeJSRegime(:,1:end-1) & supHopf(1:end-1,1:end-1); %
JSRegime = isuperlinear & esublinear & supHopf(1:end-1,1:end-1);
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
    zeroPhaseDiffregion = (reshape(hilbPhaseDiffSelected,[lI, lE]) <= 180+interval) & (reshape(hilbPhaseDiffSelected,[lI, lE]) >= 180-interval);
    
    %% JS regime plots
    f0 = figure;
    subplot(2,3,1);
    pcolor(e0_S, i0_S, 1*JSRegime); shading interp;
    colorbar;
    colormap(gca, jet); axis tight;
    subplot(2,3,2);
    pcolor(e0_S, i0_S, 1*esublinear); shading interp;
    colorbar;
    colormap(gca, jet); axis tight;
    subplot(2,3,3);
    pcolor(e0_S, i0_S, 1*isuperlinear); shading interp;
    colorbar;
    colormap(gca, jet); axis tight;
    subplot(2,3,4);
    pcolor(e0, i0, reshape(Trace,[lI,lE])); shading interp;
    colorbar;
    colormap(gca, jet); axis tight; %caxis([-0.0015, 0.0015]); 
    subplot(2,3,5);
    pcolor(e0, i0, Se); shading interp;
    colorbar;
    colormap(gca, jet); axis tight; %caxis([-0.0015, 0.0015]); 
    subplot(2,3,6);
    pcolor(e0, i0, Si); shading interp;
    colorbar;
    colormap(gca, jet); axis tight; %caxis([-0.0015, 0.0015]);
%     subplot(2,3,5);
%     pcolor(e0_S, i0_S, 1*dSe_x_e0res); shading interp;
%     colorbar;
%     colormap(gca, jet); caxis([-0.0015, 0.0015]); % axis tight;
%     subplot(2,3,6);
%     pcolor(e0_S, i0_S, 1*dSi_x_i0res); shading interp;
%     colorbar;
%     colormap(gca, jet); caxis([-0.0015, 0.0015]); % axis tight;
    %% Firing Rate plot
    
    figure('name',popln.Population_Name);
    set(gcf, 'Windowstate', 'normal');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    colormap jet
    E0 = (reshape(E0, [lI, lE])); I0 = (reshape(I0, [lI, lE]));
    subplot(121)
    pcolor(e0,i0,(reshape(eFR,[lI, lE]))); shading interp;
    hold on;
%     [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
%     contourline.LineWidth = 1.5;
%     colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('E firing rate');
    xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    
    subplot(122)
    pcolor(e0,i0,(reshape(iFR,[lI, lE]))); shading interp;
    hold on;
%     [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
%     contourline.LineWidth = 1.5;
%     colorbar;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('I Firing rate');
    xlabel('Input to E pop. (i_{E})','FontWeight','bold');
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
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('Peak Gamma Amplitude');
    
    subplot(234)
    pcolor(e0,i0,(reshape(harmonicAmp,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    title('Harmonic Amplitude');
    
    subplot(232)
    pcolor(e0,i0,(reshape(peakFreq,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('Peak Frequency (Hz)');
    
    subplot(235)
    pcolor(e0,i0,(reshape(freqRatioSelected,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    colorbar;
    axis tight;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    caxis([0 3])
    % hold on;
    % [~,contourline] = contour(e0,i0, ~isnan(reshape(freqRatioSelected,[lI, lE])), 1, 'c');
    % contourline.LineWidth = 2;
    title('Frequency Ratio (H/F)');
    
    subplot(233)
    pcolor(e0,i0,(reshape(powerRatioSelected,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    % xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    colorbar;
    axis tight;
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
    axis tight;
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
    caxis([0, ([1,1,1]*(2.^(res*[2;1;0])) ) * (2^res-1)]);
    colormap(gca, rgbmap);
%     hold on;
    % [~,contourline] = contour(e0,i0, ingoodpos & ~(zeroPhaseDiffregion), 1, 'r');
    % contourline.LineWidth = 3;
    hold off;
%     subplot(236)
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
%         cbar.TickLabels(i) = {num2str(360/(numel(cbar.TickLabels)-1)*(i-1)-180)};
%     end
    c = colorbar;
    c.Ticks = (0:6)/6; % -180:60:180; 
    c.TickLabels = {'-180','-120',' -60','  0','  60',' 120',' 180'};
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
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
    %% Gamma and Harmonic Analysis plot _ draft 8
    f2a = figure('name',[popln.Population_Name,'_1']);
    set(gcf, 'Windowstate', 'maximize');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    colormap jet
    f2b = figure('name',[popln.Population_Name,'_2']);
    set(gcf, 'Windowstate', 'maximize');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    colormap jet
    
    figure(f2a);
    subplot(231)
    pcolor(e0,i0,(reshape(peakAmp,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    caxis([0 70e-3]);
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    title('Peak Gamma Amplitude');
    
    figure(f2b);
    subplot(231)
    pcolor(e0,i0,(reshape(harmonicAmp,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    caxis([0 13e-3]);
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
    title('Harmonic Amplitude');
    
    figure(f2a);
    subplot(232)
    pcolor(e0,i0,(reshape(peakFreq,[lI, lE]))); shading interp;
    hold on;
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    colorbar;
    axis tight;
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    title('Peak Frequency (Hz)');
    
%     subplot(235)
%     pcolor(e0,i0,(reshape(freqRatioSelected,[lI, lE]))); shading interp;
%     hold on;
%     [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
%     contourline.LineWidth = 1.5;
%     set(gca, 'Color', [0.5, 0.5, 0.5]);
%     colorbar;
%     axis tight;
%     % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
%     caxis([0 3])
%     % hold on;
%     % [~,contourline] = contour(e0,i0, ~isnan(reshape(freqRatioSelected,[lI, lE])), 1, 'c');
%     % contourline.LineWidth = 2;
%     title('Frequency Ratio (H/F)');
%     
%     subplot(233)
%     pcolor(e0,i0,(reshape(powerRatioSelected,[lI, lE]))); shading interp;
%     hold on;
%     [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
%     contourline.LineWidth = 1.5;
%     set(gca, 'Color', [0.5, 0.5, 0.5]);
%     % xlabel('Input to E pop. (i_{E})','FontWeight','bold'); ylabel('Input to I pop. (i_{I})','FontWeight','bold');
%     colorbar;
%     axis tight;
%     % hold on;
%     % [~,contourline] = contour(e0,i0, ingoodpos, 1, 'c');
%     % contourline.LineWidth = 2;
%     title('Power Ratio (H/F)');
    %%
    figure(f2b);
    subplot(232)
    colormap(gca, hsv);
    pcolorSurface = pcolor(e0,i0,(reshape(wrapTo180(hilbPhaseDiffSelected),[lI, lE]))); % pcolorSurface.FaceColor = 'interp';
    set(gca, 'Color', [0.5, 0.5, 0.5]);
    % xlabel('iE','FontWeight','bold'); ylabel('iI','FontWeight','bold');
    colormap(gca, hsv);
    colorbar;
    axis tight;
    caxis([-180 180]);
%     %%
%     res = 5; %bits/color
%     x = 0:2^(3*res)-1; x=x(:);
%     rgbmap = [floor(x/2^(2*res))/(2^res-1), floor((x-floor(x/2^(2*res))*2^(2*res))/2^(1*res))/(2^res-1), floor((x-floor(x/2^(1*res))*2^(1*res))/2^0)/(2^res-1)];
%     hsvmap = hsv(numel(x));
%     noNanIndices = ~isnan(hilbPhaseDiffSelected)==1;
%     hsvdata = zeros(numel(hilbPhaseDiffSelected),3)/0;
%     hsvdata(noNanIndices,:) = hsvmap(round((wrapTo180(hilbPhaseDiffSelected(noNanIndices))+180)/360*(size(hsvmap,1))+1),:);
%     hexdata = ( hsvdata * (2.^(res*[2;1;0])) ) * (2^res-1);
%     % colormap(gca, rgbmap);
%     pcolorSurface = pcolor(e0, i0, reshape( hexdata(:), [lI,lE])); shading interp;
%     pcolorSurface.FaceColor = 'interp';
%     colormap(gca, rgbmap);
%     hold on;
%     % [~,contourline] = contour(e0,i0, ingoodpos & ~(zeroPhaseDiffregion), 1, 'r');
%     % contourline.LineWidth = 3;
%     hold off;
%     figure(f2b);
%     subplot(232)
%     hsvdata = hsvmap(round(hilbPhaseDiff/360*(size(hsvmap,1))+1),:);
%     hsvimg = reshape(hsvdata, [lI lE 3]);
%     scale = 5;
%     scale = 10; I = image(imresize(hsvimg, scale, 'bilinear'));
%     linearindices = transpose(round(1:(lI-1)/(lI*scale-1):lI))*ones(1,lE*scale) + ones(lI*scale,1)*round(0:(lE-1)/(lE*scale-1):lE-1)*lI;
%     noNanIndices = ~isnan(hilbPhaseDiffSelected(linearindices));
%     ax = gca; colordata = ax.Children.CData;
    hold off;
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
    [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
    contourline.LineWidth = 1.5;
    [~,contourline] = contour(e0,i0, ingoodpos & (zeroPhaseDiffregion), 1, 'k');
    contourline.LineWidth = 1.5;
    hold on;
    title('Phase Diff (H/F)');
    %% 
    if exist('nplotFig10','var')
        open([basefolder, '/', 'PhaseDifferentProxies.fig']);
        subplot(2,3,nplotFig10);
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
        [~,contourline] = contour(e0_S,i0_S, JSRegime, 1, 'w');
        contourline.LineWidth = 1.5;
        [~,contourline] = contour(e0,i0, ingoodpos & (zeroPhaseDiffregion), 1, 'k');
        contourline.LineWidth = 1.5;
        savefig(gcf, [basefolder, '/', 'PhaseDifferentProxies.fig']);
    end
    %%
    figure(f2a);
    print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'_A.png'],'-dpng','-r480');
    savefig(gcf, [plotsavedir, '/', 'GammaHarmonic_', plotsavename,'_A.fig'] ,'compact');
    figure(f2b);
    print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'_B.png'],'-dpng','-r480');
    savefig(gcf, [plotsavedir, '/', 'GammaHarmonic_', plotsavename,'_B.fig'] ,'compact');
    % print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.svg'],'-dsvg','-r480');
    % print([plotsavedir, '/', 'GammaHarmonic_', plotsavename,'.tif'],'-dtiff','-r480');
    disp(['Plots Saved to: ', what(plotsavedir).path,'\*_',plotsavename]);
    %%
    %%
    InpE_regime = E0(ingoodpos & (zeroPhaseDiffregion)); InpE_regime = InpE_regime(:);
    InpI_regime = I0(ingoodpos & (zeroPhaseDiffregion)); InpI_regime = InpI_regime(:);
    Phasediff_regime = wrapTo180(hilbPhaseDiffSelected(ingoodpos & (zeroPhaseDiffregion)));
    powerRatio_regime = reshape(powerRatioSelected,[lI, lE]); powerRatio_regime = powerRatio_regime(ingoodpos & (zeroPhaseDiffregion));
    regimevals = table(sortrows([InpE_regime, InpI_regime, Phasediff_regime, powerRatio_regime], 3)) %, 'VariableNames', {'Input to E','Input to I','PhaseDifference (H-G)', 'PowerRatio (H/G)'})
    save([Objsavedir,'/','Regimevals_',Objsavename,'.mat'], 'regimevals');
    %% plotting high power harmonic in-regime cases
    HPInpE_regime = InpE_regime(InpE_regime<3 & InpI_regime>8.5);
    HPInpI_regime = InpI_regime(InpE_regime<3 & InpI_regime>8.5);
    %% plot firing rates for input cases that fall within the regime
    if isempty(InpE_regime)
    else
        randind = randperm(length(InpE_regime), min(3,length(InpE_regime)));
        E_selectedcase = InpE_regime(randind); % Change this and execute section for visualizing other cases manually
        I_selectedcase = InpI_regime(randind);
        table([E_selectedcase(:), I_selectedcase(:)])%, 'VariableNames', {'Input to E','Input to I'})
        t = popln.EIpairs.t;
        timerange = t>=1.1 & t<=1.175;
        allregime = ingoodpos & (zeroPhaseDiffregion) & [[JSRegime;ones(1,size(JSRegime,2))<0],ones(size(JSRegime,1)+1,1)<0];
        
        newE0 = E0(allregime); newE0 = newE0(:);
        newI0 = I0(allregime); newI0 = newI0(:);
        AmpRatio = harmonicAmp./peakAmp; newAmpRatio = AmpRatio(allregime); newAmpRatio = newAmpRatio(:);
        newAmpRatiomax = sort(newAmpRatio, 'descend'); newAmpRatiomax = newAmpRatiomax(min(3,numel(newAmpRatiomax)));
%         newI0 = I0(allregime); newI0 = newI0(:);
        selectedcases = []; sel = find(newAmpRatio >= newAmpRatiomax); E_selectedcase = newE0(sel); I_selectedcase = newI0(sel);
        %%
        save([Objsavedir,'/','examples_',Objsavename,'.mat'], 'E_selectedcase','I_selectedcase')
        %%
        for i = 1:length(E_selectedcase)
            selectedcases = [selectedcases,find(E0(:)==E_selectedcase(i) & I0(:)==I_selectedcase(i))];
        end
%         selectedcases = [];
%         for i = 1:length(E_selectedcase)
%             selectedcases = [selectedcases,find(E0(:)==E_selectedcase(i) & I0(:)==I_selectedcase(i))];
%         end
        %%
        f = figure('name', 'Example Firing rates inside the estimated regime','windowstate','maximize');
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
        f = figure('name', 'All Firing rates inside the estimated regime','windowstate','maximize');
        set(f, 'DefaultLineLinewidth', 2.5);
        egregime = allregime; % & (E0<2.5)&(I0>8);
        nrows = floor(sqrt(sum(egregime,'all'))); ncols = ceil(sum(egregime,'all')/nrows);
        t = popln.EIpairs.t(timerange);
        iplot = find(egregime(:))
        size(iplot)
        for plotnum = 1:sum(egregime,'all')
            
            eR = popln.EIpairs.R(iplot(plotnum),timerange);
            iR = popln.EIpairs.R(popln.EIpairs.NNeuronGroups+iplot(plotnum),timerange);
            subplot(nrows,ncols,plotnum); plot(t, -(eR+iR)); title(num2str(iplot(plotnum))); axis tight;
        end
%         subplot(4,1,1); plot(t, cos(2*pi*30*t)+ 1/5*cos(2*pi*60*t)); title('Desired waveshape: cos(2{\pi} 30t) + 1/5 cos(2{\pi} 60t)');
%         eR = popln.EIpairs.R(selectedcases,timerange);
%         subplot(4,1,2); plot(t, eR); title('E cell firing rates')
%         iR = popln.EIpairs.R(popln.EIpairs.NNeuronGroups+(selectedcases),timerange);
%         subplot(4,1,3); plot(t, iR); title('I cell firing rates')
%         subplot(4,1,4); plot(t, eR+iR); title('E+I cell firing rates')
        xlabel('time (seconds)')
        
        
        %% save above example plots
        print([plotsavedir, '/', 'AllFiringRates_', plotsavename,'.png'],'-dpng','-r480');
        savefig(gcf,[plotsavedir, '/', 'AllFiringRates_', plotsavename, '.fig'],'compact');
        % print([plotsavedir, '/', 'ExampleFiringRates_', plotsavename,'.svg'],'-dsvg','-r480');
        % print([plotsavedir, '/', 'ExampleFiringRates_', plotsavename,'.tif'],'-dtiff','-r480');
%         %%
%         f = figure('name', 'Example Firing rates traces');
%         set(f, 'DefaultLineLinewidth', 2.5);
% %         t = t(timerange);
%         subplot(4,1,1); plot(t, cos(2*pi*30*t)+ 1/5*cos(2*pi*60*t)); title('Desired waveshape: cos(2{\pi} 30t) + 1/5 cos(2{\pi} 60t)');
%         eR = popln.EIpairs.R(selectedcases,timerange);
%         iR = popln.EIpairs.R(popln.EIpairs.NNeuronGroups+(selectedcases),timerange);
%         subplot(4,1,2); plot(t, eR(1,:)+iR(1,:)); title('E cell firing rates')
%         subplot(4,1,3); plot(t, eR(2,:)+iR(2,:)); title('I cell firing rates')
%         subplot(4,1,4); plot(t, eR(3,:)+iR(3,:)); title('E+I cell firing rates')
%         xlabel('time (seconds)')
%         
        %%
    end
    
    %% plot E firing rates in a grid for visual inspection
    if ~exist('nr','var')
        nr = 4;
    end
    if ~exist('nc','var')
        nc = 4;
    end
    figure('name','E firing rates');
    set(gcf, 'Windowstate', 'maximize');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    subplot(1,1,1);
    set(gca,'color','w');
    plot(0,0,'w');
    xlabel(['Subplots arranged by Input to E pop. (i_{E}): left-to-right =',num2str(e0(1)),' to ', num2str(e0(end))],'FontWeight','bold');
    ylabel(['Subplots arranged by Input to I pop. (i_{I}): left-to-right =',num2str(i0(1)),' to ', num2str(i0(end))],'FontWeight','bold');
    title('E+I cell firing rates');
    set(gcf, 'Defaultlinelinewidth', 2);
    
    E0ref = E0(:); I0ref = I0(:);
    plotrange = (popln.EIpairs.t>1 & popln.EIpairs.t <= 1+3*1/25); % t0 < t <= t0+n*1/f ==> n full cycles of an f Hz oscillation
    subsamplestep = 4;
    plotaxes = [];
    
    Ne = 1; Ni = 2; Nnrngrps = popln.EIpairs.NNeuronGroups;
    REplusI = popln.EIpairs.R((Ne-1)*Nnrngrps + (1:Nnrngrps), plotrange) + popln.EIpairs.R((Ni-1)*Nnrngrps + (1:Nnrngrps), plotrange);
    
    for inde = 1:subsamplestep:length(e0)% 0:nc/length(e0):nc-1/length(e0)
        for indi = 1:subsamplestep:length(i0)%0:nr/length(i0):nr-1/length(i0)
            indr = indi/length(i0)*(nr-1/length(i0)); indc = inde/length(e0)*(nc-1/length(e0));
            indr = ceil(indr); indc = ceil(indc);
            %         [indi, inde, indr, indc, (nr-indr)*nc + indc]
            ax = subplot(nr,nc, (nr-indr)*nc + indc);
            plotaxes = [plotaxes, ax];
            set(gca, 'color', [0.1, 0.1, 0.1]);
            hold on;
            plot(popln.EIpairs.t(plotrange), REplusI(E0==e0(inde)&I0==i0(indi),:));
            axis tight;
            aylim = ylim;
            ylim([-0.05, max(aylim(2), 0.5)]);
            hold off;
            xlabel('time (s)');
        end
    end
    linkaxes(plotaxes,'x');
    
    print([plotsavedir, '/', 'FiringRate grid_', plotsavename,'.png'],'-dpng','-r480');
    savefig(gcf,[plotsavedir, '/', 'FiringRate grid_', plotsavename, '.fig'],'compact');
    %% plot fft of E firing rates in a grid for visual inspection
    if ~exist('nr','var')
        nr = 5;
    end
    if ~exist('nc','var')
        nc = 5;
    end
    figure('name','FFT of E+I firing rates');
    set(gcf, 'Windowstate', 'maximize');
    set(gcf, 'color', 'w');
    set(gcf,'InvertHardcopy','off');
    subplot(1,1,1);
    set(gca,'color','w');
    plot(0,0,'w');
    xlabel(['Subplots arranged by Input to E pop. (i_{E}): left-to-right =',num2str(e0(1)),' to ', num2str(e0(end))],'FontWeight','bold');
    ylabel(['Subplots arranged by Input to I pop. (i_{I}): left-to-right =',num2str(i0(1)),' to ', num2str(i0(end))],'FontWeight','bold');
    title('FFT of E+I cell firing rates');
    set(gcf, 'Defaultlinelinewidth', 1.5);
    
    E0ref = E0(:); I0ref = I0(:);
    plotrange = (f_fft>25 & f_fft <= 150); % t0 < t <= t0+n*1/f ==> n full cycles of an f Hz oscillation
    subsamplestep = 4;
    plotaxes = [];
    for inde = 1:subsamplestep:length(e0)% 0:nc/length(e0):nc-1/length(e0)
        for indi = 1:subsamplestep:length(i0)%0:nr/length(i0):nr-1/length(i0)
            indr = indi/length(i0)*(nr-1/length(i0)); indc = inde/length(e0)*(nc-1/length(e0));
            indr = ceil(indr); indc = ceil(indc);
            %         [indi, inde, indr, indc, (nr-indr)*nc + indc]
            ax = subplot(nr,nc, (nr-indr)*nc + indc);
            plotaxes = [plotaxes, ax];
            set(gca, 'color', [0.1, 0.1, 0.1]);
            hold on;
            plot(f_fft(plotrange), abs(-R_fft(E0==e0(inde)&I0==i0(indi), plotrange)));
            axis tight;
            aylim = ylim;
            ylim([-0.5, max(aylim(2), 1)]);
            hold off;
            xlabel('frequency (Hz)');
        end
    end
    linkaxes(plotaxes, 'x');
    
    print([plotsavedir, '/', 'FFT grid_', plotsavename,'.png'],'-dpng','-r480');
    savefig(gcf,[plotsavedir, '/', 'FFT grid_', plotsavename, '.fig'],'compact');
    %%
    Population_Name = popln.Population_Name;
    R = popln.EIpairs.R;
    t = popln.EIpairs.t;
    W = popln.EIpairs.W;
    fftLFP = R_fft;
    fftf = f_fft;
    tau = popln.EIpairs.tau;
    Input = popln.EIpairs.Input;
    DynamicInputParameters = popln.DynamicInputParameters;
    save([Objsavedir,'/','PoplnObject_',Objsavename,'.mat'], 'Population_Name','R', 't','fftLFP','fftf', 'W', 'tau', 'Input', 'DynamicInputParameters','-v7.3');
    
end