function [graphInfo] = plotGammaHarmonicSummary(graphInfo, SummaryPlotDim, InRegime, rowparametersDescription, rowvalues, colparametersDescription, colvalues)

    %%
    assert(exist('rowparametersDescription','var')==1, 'Provide simparameterDescription');
    %%
    if (~exist('simparametersDescription','var'))
        simparametersDescription = ['<row_parameter_symbol>']
    end
    if (~exist('rowparametersDescription','var'))
        rowparametersDescription = simparametersDescription;
    end

    if (~exist('colparametersDescription','var'))
        colparametersDescription = ''; %'<col_parameter_symbol>';
    end
    %%
    assert(~isempty(SummaryPlotDim),'Specify number of rows to be plotted in SummaryPlotDim');    % simparameterlist is a cell containing arrays of min and max parameter values
    RowSizeW = 10; RowSizeH = (SummaryPlotDim==2)*RowSizeW + (SummaryPlotDim==1)*1;% cm
    WidHt = [RowSizeW,RowSizeH];
    figbyaxScale = [1.2, 1.2]; 
    %%
    if (isempty(graphInfo))

        SummFig = figure('name','Summary Figure', 'WindowState','Normal');
        %%
        figWidHt = figbyaxScale.*WidHt;
        prevFigHt = 0;
    else
        SummFig = graphInfo.fig;
        prevFigHt = graphInfo.prevFigHt;
        figWidHt = figbyaxScale.*WidHt + [0,prevFigHt];
    end
    set(SummFig,'Resize','off');
    set(SummFig, 'PaperUnits', 'centimeters');
    set(SummFig, 'PaperSize', figWidHt);
    set(SummFig, 'PaperPositionMode', 'manual');
    set(SummFig, 'PaperPosition', [[0 0], figWidHt]);
    set(SummFig,'units','centimeters');
%     set(SummFig,'position',[[0,0],figWidHt]);
    axposition = [[0.75, 0.9].*(1-1./figbyaxScale) + [0,prevFigHt/figWidHt(2)], [1,1-prevFigHt/figWidHt(2)]./figbyaxScale];
    figscale = prevFigHt/figWidHt(2)*[1 1 1 1];
    for ax = SummFig.Children
        set(ax, 'units', 'Normalized');
        set(ax, 'fontunits', 'Normalized');
    end
    if prevFigHt > 0
        set(SummFig, 'Position', SummFig.Position.*figscale);
    end
    for ax = SummFig.Children
        set(ax, 'units', 'centimeters');
        set(ax, 'fontunits', 'centimeters');
    end
    if prevFigHt > 0
        set(SummFig, 'Position', SummFig.Position.*[1,1,1,1/figscale(end)]);
    end
    subplot('Position',axposition);
    set(gca,'units','centimeters');
    camva('manual'); daspect manual; pbaspect manual;
    
    assert(exist('rowvalues','var')==1, 'The values of the 1^st varied parameter are required to form the x-axis of plot')
    if SummaryPlotDim == 1
        yticks(gca,[]);
        colvalues = 0.5*ones(size(rowvalues));
    else
        assert(exist('colvalues','var')==1, 'The values of the 2^nd varied parameter are required to form the y-axis of plot')
    end

    if (~exist('TrueColor','var'))
        TrueColor = [0, 0.95, 0];
    end
    if (~exist('FalseColor','var'))
        FalseColor = [0.66, 0.33, 0.33];
    end
    colors = InRegime(:)*TrueColor + (1-InRegime(:))*FalseColor;
    scatter(rowvalues,colvalues,[], colors,'filled');

    axis tight;
    rowlims = [min(rowvalues), max(rowvalues)];
    try
        xlim(rowlims + [-1,+1]*abs(diff(rowlims))/7.5*0.5);
    catch
    end
    if SummaryPlotDim == 1
        ylim([0, 1]);
    else
        collims = [min(colvalues), max(colvalues)];
        try
            ylim(collims + [-1,+1]*abs(diff(collims))/7.5*0.5);
        catch
        end
    end
    if SummaryPlotDim == 1
        yticks(gca,[]);
    end
    xlabel(rowparametersDescription,'FontWeight','bold');
    ylabel(colparametersDescription,'FontWeight','bold');
    %%
    graphInfo.fig = SummFig;
    graphInfo.prevFigHt = SummFig.PaperPosition(end);
    SummFig.Position
%     set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
end