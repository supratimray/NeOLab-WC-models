classdef GenerateNeuronPopulation < handle & dynamicprops & matlab.mixin.Copyable
    properties
        % Properties could be added to an EIpopulation instance after creation using 'addprops' function.
        % Functions for assigning weights and Inputs can be created outside the class and matched to the handles Wfunction and SetEIInput handles respectively
        
        EIpairs = [], ids = [1], Ngroups = 1; NStateVars = 2;
        %% User variables for Population 
        Dim = [1], Population_Name='';
        %% Connectivity parameters and Functions
        circular=[]; %Mark the dimension along which the pairs are connected in a ring (say, for modelling Normalization)
        connectivitySigma = 0;
        connectivityMaxW = [0.1, 0.089 ; 0.38, 0.096];
        Wfunction = @(obj, nvars) obj.GaussianWtDist(obj.connectivitySigma, obj.connectivityMaxW, nvars);
        %% Input, Initial state and Simulation details
%         SetEIInput = @(obj, tspan, Input) obj.input(tspan,Input);
        R0 = 0;
        tspan = [0, 2500]*1e-3;
        EI_update;
        % For DYNAMIC NEURONAL INPUT/PROPERTY conditions in simulations
        DynamicInputParameters = []; % set/subset of neuronal inputs/ parameters controlling inputs that are updated during ode simulation
        NDynamicInputParameters = 0 % Default: 0-static stimulus condition
        % A separate function must be written to compute Inputs to neurons from NDynamicInputParameters. 
        % The aforementioned @(DynParameters)->NeuronalInputs function may be used within a InputUpdatefn, 
                % which will be passed as an argument to obj.EI_update_Dynamicinput (method defined below) 
                % if obj.EI_update = @(t,y)EI_update_Dynamicinput(t,y, @(t,in_latest)InputUpdatefn(t,in_latest))
        UpdateInputs = @(obj, t) []; % Modify function to Instantaneously Vary inputs. Must set input values for each instant and return the same, for retrospective tracing/plotting (after simulation)
    end
    methods
        function setupEIPopulation(obj, dimensions, NStateVars, name, tau, Wfunction, setSummationNonlinearity)
            obj.Dim = dimensions;
            
            if ~exist('NStateVars', 'var')
                NStateVars = 2; % Defaults to EIpair network
                disp('Setting 2 state variables per Neuron group.')
            end
			obj.NStateVars = NStateVars;
            
            obj.Ngroups = 1;
            for i = 1:length(dimensions)
                obj.Ngroups = obj.Ngroups *dimensions(i);
            end
            
            if exist('name', 'var')
                obj.Population_Name = name;
            end
            if exist('Wfunction','var')
                obj.Wfunction = Wfunction;
            end
            if exist('setSummationNonlinearity','var')
                setSummationNonlinearity(obj.EIpairs); % function handle @(obj.EIpairs) -> adds required attributes to the argument class, and sets obj.EIpairs.InSum
            end
%             disp(obj.ids)
            if isrow(obj.connectivitySigma)
                obj.connectivitySigma = transpose(obj.connectivitySigma);
            end
            obj.EIpairs = [];
            obj.ids = [];
            if iscolumn(dimensions)
                dimensions = transpose(dimensions);
            end
            obj.ids = transpose(0:obj.Ngroups-1)*ones(1,length(dimensions));
            cumdims = cumprod(dimensions);
            tempid = (transpose(1:obj.Ngroups)*ones(1,size(dimensions,2)) -1) ;
            tempid = floor(tempid / diag([1,cumdims(1:end-1)]))   ;         
            obj.ids = 1 + round(tempid - floor(tempid/ diag(dimensions))* diag(dimensions));
            if ~ exist('tau','var')
                obj.EIpairs = Neurons(dimensions, NStateVars, @()obj.Wfunction(obj, NStateVars));
            else
                obj.EIpairs = Neurons(dimensions, NStateVars, @()obj.Wfunction(obj, NStateVars), tau);
            end
            obj.EI_update = @(t,y) obj.EI_update_StaticInput(t,y);
        end
                
        
        function yval = findIndex(obj, y,i) % TODO push inside EIpairarray class with Index dist function
            yval = [y(i, :); y(i+obj.Ngroups,:)];
        end
        
        function input(obj,Is, tspan, R0) % TODO keep this here
            if exist('tspan','var')>0
                obj.tspan = tspan;
            end
                    
            if exist('Is','var')
                Is = Is(:);
                if ~isempty(Is)
                    if numel(obj.tspan)==2
                        ti = obj.tspan(1); tf = obj.tspan(2);
                        obj.tspan = ti : min(1,tf-ti)*1e-4 : tf;
                    end
                    obj.EIpairs.tspan = obj.tspan;

                    obj.EIpairs.Input = zeros(obj.Ngroups*obj.EIpairs.NStateVarsPerGroup,1);
                    ln = obj.Ngroups;
                    for i = 0:obj.EIpairs.NStateVarsPerGroup-1
                        obj.EIpairs.Input(i*ln+1:(i+1)*ln) = Is(ceil((i*end+1)/obj.EIpairs.NStateVarsPerGroup):ceil((i+1)*end/obj.EIpairs.NStateVarsPerGroup));
                    end
                    
                end
            end
            if exist('R0','var')
                obj.R0 = R0;
            end
        end
        
         function wts = EIGaussianWtDist(obj, sigma, K, ~)
             delxsqEE = zeros(obj.Ngroups);
             delxsqEI = zeros(obj.Ngroups);
             delxsqIE = zeros(obj.Ngroups);
             delxsqII = zeros(obj.Ngroups);
             for coor = 1:length(obj.Dim)
                 temp = obj.ids(:,coor)*ones(1,obj.Ngroups);
                 normdiff = (temp - transpose(temp));
                 if obj.circular(coor) == 1
                     normdiff2 = obj.Dim(coor) - abs(normdiff);
                     normdiff = min(abs(normdiff), normdiff2);
                 end
                 normdiff = normdiff.^2;
                 delxsqEE = delxsqEE + normdiff/(sigma(coor, 1,1)^2);
                 delxsqEI = delxsqEI + normdiff/(sigma(coor, 1,end)^2);
                 delxsqIE = delxsqIE + normdiff/(sigma(coor, end,1)^2);
                 delxsqII = delxsqII + normdiff/(sigma(coor, end,end)^2);
             end
             wtE = [K(1,1)*exp(-0.5*delxsqEE),K(1,2)*exp(-0.5*delxsqEI)];
             wtI = [K(2,1)*exp(-0.5*delxsqIE),K(2,2)*exp(-0.5*delxsqII)];
             wts = [wtE; wtI];
         end
         
         function wts = GaussianWtDist(obj, sigma, K, nStatevars) % nStatevars = 2 for EI pairs, 3 for EIS units and so on ==> # state variables per neuron group/column
             delxsq = cell(nStatevars, nStatevars);
             dzer = ones(obj.Ngroups);
             for i = 1:nStatevars
                for j = 1:nStatevars
                    delxsq{i,j} = dzer;
                end
             end
             for coor = 1:length(obj.Dim)
                 temp = obj.ids(:,coor)*ones(1,obj.Ngroups);
                 normdiff = (temp - transpose(temp));
                 if ~isempty(obj.circular)
                     if obj.circular(coor) == 1
                        normdiff2 = obj.Dim(coor) - abs(normdiff);
                        normdiff = min(abs(normdiff), normdiff2);
                     end
                 end
                 % expnormdiff = exp(-normdiff.^2);
                 for i = 1:nStatevars
                    for j = 1:nStatevars 
                        xyz = ceil([coor/size(sigma,1), i/size(sigma,2), j/size(sigma,3)]/nStatevars);
                        delxsq{i,j} = delxsq{i,j} .* (exp(-normdiff.^2).^(1/sigma(xyz(1), xyz(2), xyz(3))^2)); % Change this portion for a nonlinearity/stochasticity
                    end
                 end
             end
             wts = [];
             for i = 1:nStatevars
                wt = [];
                for j = 1:nStatevars
                    wt = [wt, K(i,j)*delxsq{i,j}]; % Change this portion for a nonlinearity/stochasticity
                end
                wts = [wts; wt];
             end
         end
         
        function wt = weightDistribution(obj, pair_n, sigma, K)
            % SSN model
            if exist('K','var') == 0
                K = [1.05,0.5;0.5,0.1]; %%% Use SSN maxweights
            end
            if exist('sigma', 'var') == 0
                sigma = 0;
            end
            
            wt = zeros(2, 2*obj.Ngroups);
            wt(:, pair_n + [0, obj.Ngroups]) = K;
            
            if any(sigma > 0)
                for i = 1:obj.Ngroups
                    if i ~= pair_n
                        %%%
                        x1 = obj.ids(pair_n,:); x2 = obj.ids(i,:);
                        %exp( -0.5 * sum(((x1-x2)./sigma).^2))
                        delx = x1-x2; delxp = obj.Dim - abs(delx);
%                         abs([delx;delxp])
%                         abs([delx(obj.circular == 1);delxp(obj.circular == 1)])
                        delx(obj.circular == 1) = min(abs(delx(obj.circular==1)), delxp(obj.circular==1));
%                         [delx(obj.circular == 1);abs(delxp(obj.circular == 1))]
                        for ic = 0:1
                            for jc= 0:1
                                wt(ic+1, i+jc*obj.Ngroups) = K(ic+1,jc+1) *exp( -0.5 * sum((transpose(delx)./sigma(:,(1-ic)+ic*size(sigma,2),(1-jc)+jc*size(sigma,3))).^2));
%                                 ic, jc, sigma(:,(1-ic)+ic*size(sigma,2),(1-jc)+jc*size(sigma,3))
                            end
                        end
                    end
                end
            end
        end
        
        
		function ode(obj, solver, odeopts)            
			overwriteFlag = isempty(obj.EIpairs.t) || obj.tspan(1) < obj.EIpairs.t(end);
            if overwriteFlag
                obj.EIpairs.R0 = zeros(obj.EIpairs.NNeuronGroups*obj.EIpairs.NStateVarsPerGroup,1);
                for i = 0:obj.EIpairs.NStateVarsPerGroup-1
                    obj.EIpairs.R0(i*obj.EIpairs.NNeuronGroups+1:(i+1)*obj.EIpairs.NNeuronGroups) = obj.R0(ceil((i*end+1)/obj.EIpairs.NStateVarsPerGroup):ceil((i+1)*end/obj.EIpairs.NStateVarsPerGroup));
                end
            else
                obj.EIpairs.R0 = obj.EIpairs.R(:,end);
            end
            
            disp('Simulation begun')   
            if ~exist('odeopts','var')         
                odeopts = odeset('RelTol',1e-3,'AbsTol',1e-3);
            end
            if ~exist('solver','var')
                solver = [];
            end
            if isempty(solver)
                disp('Using native ode solver');
                [ts, ys] = ode113(@(t,y) obj.EI_update(t,y), obj.tspan, [obj.EIpairs.R0; obj.DynamicInputParameters], odeopts);
            else
                disp('Using specified solver');
                [ts, ys] = solver(@(t,y) obj.EI_update(t,y), obj.tspan, [obj.EIpairs.R0; obj.DynamicInputParameters]);
            end
            t = [];
            if overwriteFlag
                t = transpose(ts);
                obj.EIpairs.t = []; obj.EIpairs.R = [];
            else
                t = [t, transpose(ts)];
            end
            
            tau_is_0 = (obj.EIpairs.tau==0)==1;
            obj.EIpairs.t = [obj.EIpairs.t, t];
            ys_tau0corrected = transpose(ys(:,1:end-obj.NDynamicInputParameters));
            
            if isempty(obj.UpdateInputs(obj,t(1)))
                inputs = obj.EIpairs.Input;
            else
                inputs = [];
            end
            
            if (any(tau_is_0)) || ~isempty(obj.UpdateInputs(obj,t(1)))
                
                for col = 1:length(t)
                    Is = obj.UpdateInputs(obj, t(col));
                    inputs = [inputs, Is(:)];
                    
                    obj.input(Is);
                    summedVal = obj.EIpairs.InSum(...
                                obj.EIpairs, ...
                                t(col), ...
                                ys_tau0corrected(:,col), ...
                                nan ... % recordflag = nan -> signals external Insum function not to record these inputs/noises (or even use pre-recorded input values) since these are revision routines to reconstruct instantaneous variables
                                );
                    ys_tau0corrected(tau_is_0,col) = summedVal(tau_is_0);
                end
            end
            obj.EIpairs.R = [obj.EIpairs.R, ys_tau0corrected];
            obj.DynamicInputParameters = transpose(ys(:,end-obj.NDynamicInputParameters+1:end));
            obj.EIpairs.Input = inputs;
            disp('Simulation complete!')
%             for i = 1:obj.Ngroups
% 				if overwriteFlag
% 					obj.t = transpose(ts);
% 					obj.EIpairs(i).t = []; obj.EIpairs(i).R = [];
% 				else
% 					obj.t = [obj.t, transpose(ts)];
% 				end
% 				
% 				obj.EIpairs(i).t = [obj.EIpairs(i).t, obj.t];
% 				obj.EIpairs(i).R = [obj.EIpairs(i).R, transpose(ys(:, i + [0,obj.Ngroups]))];
% 			end
        end
        
        function dR_dt = EI_update_StaticInput(obj, t, y)  % 'y' is a [2*#EIpairs]x1 column vector of activities [<Exc1;2;;;>; <Inh1;2;;;>]
            objArr = obj.EIpairs;
            obj.input(obj.UpdateInputs(obj, t));
            y_latest = y(:,end);
%             size(y_latest),
%             y_latest(1), 
%             t,
            tau_is_0 = (objArr.tau==0);
            y_latest_tau0corrected = y_latest;
            recordflag = 1; % incase the external summation function intends to record the proper input values (useful when noise is added at the summation function or the update input function)
            summedVal = objArr.InSum(objArr, ...
                        t(end), ...
                        y_latest_tau0corrected , recordflag); 
            y_latest_tau0corrected(tau_is_0) = summedVal(tau_is_0); % updates instantaneous variables (recordflag = 1)
            recordflag = 2;
            summedVal = objArr.InSum(objArr, ...
                        t(end), ...
                        y_latest_tau0corrected , recordflag); % updates the non-instantaneous variables (recordflag = 2)
            dR_dt = zeros(size(y_latest_tau0corrected));
            dR_dt(~tau_is_0) = (-y_latest(~tau_is_0) ...
                                +summedVal(~tau_is_0) ...
                                ) ./objArr.tau(~tau_is_0);
            dR_dt(tau_is_0) = 0;
            t_intspan = any(obj.tspan==t) | any(abs([obj.tspan(1),obj.tspan]+[obj.tspan,obj.tspan(end)] - 2*t) <= 0.00002);
%             t_del = (round(t*1000) - floor(t*10)*100);
            if t_intspan
                fprintf("\tsimulation @ %0.3f s\n",t)
            end
%             disp([dR_dt, Rlatest, [dR_dt(1).*objArr(:).tau(1);dR_dt(2).*objArr(:).tau(end)]+Rlatest])
        end
        
        %% TOCHECK: TODO: modify function to provide dynamic input in experiments
        function dy_dt = EI_update_DynamicInput(obj, t, y, InputUpdatefn) 
            % 'y' is a [2*#EIpairs+#Inputs]x1 column vector of activities [<Exc1;2;;;>; <Inh1;2;;;>]
            % 'inputUpdatefn' is a function handle @(t, in_latest) that- Sets objArr.Input property as per in_latest & Calculates/returns the next update/sample of the Input 
            objArr = obj.EIpairs;
            obj.input(obj.UpdateInputs(t));
            ninputs = obj.NDynamicInputParameters;

            y_latest = y(1:end-ninputs,end);
            in_latest = y(end-ninputs+1:end,end);
%             size(y_latest),
%             y_latest(1), 
%             t,
            dIn_dt = InputUpdatefn(t, in_latest);
            tau_is_0 = (objArr.tau==0);
            recordflag = 1;
            y_latest_tau0corrected = y_latest;
            summedVal = objArr.InSum(objArr, ...
                        t(end), ...
                        y_latest_tau0corrected , recordflag); % updates instantaneous variables
            y_latest_tau0corrected(tau_is_0) = summedVal(tau_is_0);
            recordflag = 2;
            summedVal = objArr.InSum(objArr, ...
                        t(end), ...
                        y_latest_tau0corrected , recordflag); % updates non-instantaneous variables 
            dR_dt = zeros(size(y_latest_tau0corrected));
            dR_dt(~tau_is_0) = (-y_latest(~tau_is_0) ...
                                +summedVal(~tau_is_0) ...
                                ) ./objArr.tau(~tau_is_0);
            dR_dt(tau_is_0) = 0; % derive in calcR_tau0 function
            dy_dt = [dR_dt; dIn_dt];
%             disp([dR_dt, Rlatest, [dR_dt(1).*objArr(:).tau(1);dR_dt(2).*objArr(:).tau(end)]+Rlatest])
            t_del = (round(t*1000) - floor(t*100)*10);
            if t_del == 9 | t_del ==0
                disp(t)
            end
        end

        %% Create time frequency analysis of activity of types, 'Power', 'Magnitude', 'Phase'
		function [f,Ets,Its] = timeFreq(type, index)
				
		end
		
		%% Modify nullclines functions
		
    end
end