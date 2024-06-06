function [NC,Quiver] = describeDynamics(fig, ax, model, modelinputs, NCvariableIDs, variableBoundsXY, nPts, ncColorXY, nQuiverPts, quiverColorXY)

    if ~exist('NCvariableIDs','var')
        NCvariableIDs = 1:2;
    end
    if ~exist('nPts','var')
        nPts = [];
    end
    if isempty(nPts)
        nPts = 1e2;
    end    
    if ~exist('nQuiverPts','var')
        nQuiverPts = [];
    end
    if isempty(nQuiverPts)
        nQuiverPts = 15;
    end    
    
    if ~exist('ncColorXY','var')
        ncColorXY = {'b','r'};
    end
    if ~exist('quiverColorXY','var')
        quiverColorXY = 'g';
    end

    %% Quiver plots
    nstatevars = sum(NCvariableIDs>0);
    
    plotvars = (1:nstatevars); 
    plotvars = plotvars(NCvariableIDs);
    
    % Stage 1 : candidate state values for finding nullclines
    vargrid = cell(nstatevars,1);
    step = 1e100;
    for i=plotvars(:)'
        interval = variableBoundsXY{i};
        interval = diff([interval(1), interval(end)]);

        if interval>0
            step = min(step,min(interval)/nQuiverPts);
        else
            error(['Invalid value bounds specified for State variable #', i, 'of model. Only single value specified']);
        end
    end
    for i=1:nstatevars
        vargrid{i} = linspace(variableBoundsXY{i}(1),variableBoundsXY{i}(2), nQuiverPts);
        if ~any(plotvars==i)
            vargrid{i} = variableBoundsXY{i}(1);
        end
    end
    
    [varCombinations, ~, nsims] = unwrapParameters(vargrid);
    varVals = varCombinations(:);
    
    
    unwrappedinput = repmat(modelinputs(:)',nsims,1);

    modelpop = model(nsims);
    modelpop.input(unwrappedinput(:), 0:1e-3:2);
    
%     rhs = @(y) (-y + modelpop.EIpairs.InSum(modelpop.EIpairs,0, y, false))./modelpop.EIpairs.tau(:);
    rhs = @(y) modelpop.EI_update(1000,y(:));
    dr_dt = reshape( rhs(varVals), [], nstatevars);
    dr_dt_XY = dr_dt(:, plotvars(:)');
    i_XY     = varCombinations(:,plotvars(:)');
    
    
    figure(fig); subplot(ax);
    hold on;
    quiver(i_XY(:,1), i_XY(:,2), dr_dt_XY(:,1), dr_dt_XY(:,2),'color',quiverColorXY);
    hold on;
    Quiver.XY = i_XY; Quiver.UV = dr_dt_XY;

    %% NullClines
    % variableStatesIDs can be logical or integer indices
    nstatevars = sum(NCvariableIDs>0);
    
    plotvars = (1:nstatevars); 
    plotvars = plotvars(NCvariableIDs);
    
    % Stage 1 : candidate state values for finding nullclines
    vargrid = cell(nstatevars,1);
    step = 1e100;
    for i=plotvars(:)'
        interval = variableBoundsXY{i};
        interval = diff([interval(1), interval(end)]);

        if interval>0
            step = min(step,min(interval)/nPts);
        else
            error(['Invalid value bounds specified for State variable #', i, 'of model. Only single value specified']);
        end
    end
    for i=1:nstatevars
        vargrid{i} = linspace(variableBoundsXY{i}(1),variableBoundsXY{i}(2), nPts);
        if ~any(plotvars==i)
            vargrid{i} = variableBoundsXY{i}(1);
        end
    end
    
    [varCombinations, ~, nsims] = unwrapParameters(vargrid);
    varVals = varCombinations(:);
    
    
    unwrappedinput = repmat(modelinputs(:)',nsims,1);

    modelpop = model(nsims);
    modelpop.input(unwrappedinput(:), 0:1e-3:2);
    
    residual = @(y) modelpop.EI_update(10000,y(:));% @(y) (-y + modelpop.EIpairs.InSum(modelpop.EIpairs,0, y, false))./modelpop.EIpairs.tau(:);
    
    rmserror = reshape((residual(varVals).^2).^0.5,[],nstatevars);
    
    minplotvars = min(plotvars,[],"all");
    maxplotvars = max(plotvars,[],"all");
    
    rmserror_min = reshape(rmserror(:,minplotvars), nPts, nPts); 
    focusid_min = (rmserror_min <= prctile(rmserror_min,100*1.5/nPts,1)); % Find complementary state variable that causes current variable's nullcline to vanish for each value of current variable
    rmserror_max = reshape(rmserror(:,maxplotvars), nPts, nPts); 
    focusid_max = (rmserror_max <= prctile(rmserror_max,100*1.5/nPts,2));
    
    ncxy = cell(1,nstatevars);
    ncxy{minplotvars} = varCombinations(focusid_min(:), :);
    ncxy{maxplotvars} = varCombinations(focusid_max(:), :);
    
    Xno = plotvars(1); Yno = plotvars(2);
    ncxy_X = ncxy{Xno};
    ncxy_Y = ncxy{Yno};
    varVals_XYid = [zeros(size(ncxy_X))+Xno; zeros(size(ncxy_Y))+Yno];    
    Xid = varVals_XYid(:,1)==Xno; Yid = varVals_XYid(:,1)==Yno;
%     varVals_stage1result = [ncxy_X; ncxy_Y];
    
    % Stage 2 : Converging on nearest nullcline point
    step = step/100;
    
    nsims = size(varVals_XYid,1);
    unwrappedinput = repmat(modelinputs(:)',nsims,1);

    modelpop = model(nsims);
    modelpop.input(unwrappedinput(:), 0:1e-3:2);
    
    residual = @(y) modelpop.EI_update(10000,y(:));% @(y) (-y + modelpop.EIpairs.InSum(modelpop.EIpairs,0, y, false))./modelpop.EIpairs.tau(:);
    
    for i=1:10000
        dt = 1/1000;
        
        varVals = [ncxy_X; ncxy_Y];
        y = varVals(:);
        rhs = reshape(residual(y),[],nstatevars);

        % computing partial derivatives
        stepVals = zeros(size(varVals,1),nstatevars);
        stepVals(:,Xno) = step*Yid; % Perturb y in x-nullcline equations
        stepVals(:,Yno) = step*Xid; % Perturb x in y-nullcline equations
        dely = stepVals.*randn(size(stepVals));
        dely(dely==0) = stepVals(dely==0);
        gradient = reshape(residual(y(:)+dely(:))-residual(y(:)),[],nstatevars);

        dy_update = 0*varVals;
        % iterative convergence to nullcline: d((dy/dt)^2)/dt = - (dy/dt)^2
        % For a given nullcline variable 'compns', only its complementary variable 'ns' 
        % is updated at each iteration
        for ns=1:nstatevars
            if any(plotvars==ns)
                compns = plotvars(plotvars~=ns);
                nc_compns = varVals_XYid(:,1) == compns;
                % computing partial derivative of "complementary variable's
                % ('compns')'s time-derivative" in the complementary variable ('compns') nullcline
                % w.r.t variable corresponding to 'ns' 
                dnccompns_ns = gradient(nc_compns,compns)./dely(nc_compns,ns);
                dnccompns_ns(dnccompns_ns.^2<1e-100) = 1e-10;  
                
                dy_update(nc_compns, ns) = dt * (-0.5*rhs(nc_compns,compns)./dnccompns_ns);
            end
        end
        ncxy_X = ncxy_X + dy_update(Xid,:);
        ncxy_Y = ncxy_Y + dy_update(Yid,:);
    end

    NC.Xnullcline=ncxy_X(:,plotvars);
    NC.Ynullcline=ncxy_Y(:,plotvars);
    % plot nullclines;
    figure(fig); subplot(ax);
    hold on;
    scatter(ncxy_X(:,plotvars(1)),ncxy_X(:,plotvars(2)),[],ncColorXY{1},'marker','.');
    hold on;
    scatter(ncxy_Y(:,plotvars(1)),ncxy_Y(:,plotvars(2)),[],ncColorXY{2},'marker','.');
    xlim(variableBoundsXY{plotvars(1)});
    ylim(variableBoundsXY{plotvars(2)});
end