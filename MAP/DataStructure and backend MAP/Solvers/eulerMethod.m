function [tVals,yVals] = eulerMethod(eqnName,tVals,y0,wcParams,stimParams) %#ok<INUSD>

dt = tVals(2)-tVals(1);
yVals(:,1)=y0;
for i=2:length(tVals)
    y=yVals(:,i-1);
    if exist('wcParams','var') && exist('stimParams','var')
        dy = eval([eqnName '([],y,wcParams,stimParams);']);
    else
        dy = eqnName(tVals(i), y);
    end
yVals(:,i) = yVals(:,i-1)+dy*dt;
end
tVals=tVals';
yVals=yVals';