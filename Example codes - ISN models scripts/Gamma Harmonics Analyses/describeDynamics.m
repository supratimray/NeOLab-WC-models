% Given an equation, generate a phase-plane diagram along with nullclines
% and fixed points

function describeDynamics(eqnName,xList,yList,wcParams,stimParams) %#ok<INUSD>

% Generate Quiver Plot
[x,y] = meshgrid(xList,yList);
u = zeros(size(x));
v = zeros(size(x));
eNC = zeros(size(x));
iNC = zeros(size(x));
fixedPoints = zeros(numel(x),2);
options = optimset('Display','off','MaxFunEvals',500); %#ok<*NASGU>
for i=1:numel(x)
    dy = eval([eqnName '([],[x(i) y(i)],wcParams,stimParams);']);
    u(i) = dy(1); v(i) = dy(2);
    
    eNC(i) = eval(['fzero(@(I)' eqnName '(0,[x(i) I],wcParams,stimParams,1),y(i),options)']);
%    eNC(i) = eval(['fzero(@(E)' eqnName '(0,[E y(i)],wcParams,stimParams,1),x(i),options)']); % This works equally well
    iNC(i) = eval(['fzero(@(I)' eqnName '(0,[x(i) I],wcParams,stimParams,2),y(i),options)']);
    fixedPoints(i,:) = eval(['fsolve(@(EI)' eqnName '(0,EI,wcParams,stimParams),[x(i) y(i)],options)']);
end

uniqueFixedPoints = uniquetol(fixedPoints,1e-4,'ByRows',true);

% plot nullclines.
xCol=x(:); eNC=eNC(:); iNC=iNC(:);
[xSort,I] = sort(xCol);
eNCSort = eNC(I); iNCSort = iNC(I);

plot(xSort,eNCSort,'ro'); hold on;
plot(xSort,iNCSort,'bo');

% Trajectory from origin
tVals=1:1000;
y0 = [0 0];
[~,trajectory] = eval(['ode45(@(t,y)' eqnName '(t,y,wcParams,stimParams),tVals,y0)']);
plot(trajectory(:,1),trajectory(:,2),'k');
legend('dE/dt=0','dI/dt=0','trajectory');

% Quiver plot
quiver(x,y,u,v,'g');
xlabel('E'); ylabel('I'); 
axis([min(xList) max(xList) min(yList) max(yList)]);


% Fixed points
numFixedPoints  = size(uniqueFixedPoints,1);

for i=1:numFixedPoints
    [~,~,~,~,Jacobian] = eval(['fsolve(@(EI)' eqnName '(0,EI,wcParams,stimParams),uniqueFixedPoints(i,:),options)']);
    eigs = eig(Jacobian);
    L1 = eigs(1); L2=eigs(2);
    
    if (L1==0 || L1==0)
        % Work on these degenerate cases
    elseif (isreal(L1) && isreal(L2))
        if max(L1,L2)<0
            textLabel = 'stable node';
        elseif min(L1,L2)>0
            textLabel = 'unstable node';
        else
            textLabel = 'saddle';
        end
    else
        if real(L1)<0
            textLabel = 'stable focus';
        elseif real(L1)>0
            textLabel = 'unstable focus';
        else
            textLabel = 'center';
        end
    end
    
    plot(uniqueFixedPoints(i,1),uniqueFixedPoints(i,2),'ko'); hold on;
    text(uniqueFixedPoints(i,1)+0.01,uniqueFixedPoints(i,2)+0.01,textLabel);
end
end