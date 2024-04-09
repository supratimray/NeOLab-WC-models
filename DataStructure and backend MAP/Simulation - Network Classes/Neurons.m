% Excitatory-Inhibitory population generation and simulations

classdef Neurons < handle & dynamicprops
    properties
        ids = [1];
        NNeuronGroups = 1; 
        % No of neuron groups (e.g., EI-pairs if NStateVarsPerGroup is 2)
        NStateVarsPerGroup = 2; 
        % No of state variables (updated over time, e.g, E firing rate, I firing rate, Adaptation, Surround suppression) 
        % per NeuronGroup (default: 2- e.g. EI pair)
        Dim = [1];
        R; dR_dt; R0 = [0;0]; tau=[20e-3; 10e-3]; Input = 20;
		W = [1 -1; 1 -1]; % Generally a 2x{#EIpairs} matrix containing the absolute weights of [<Excitatory 1...n>, <Inhibitory 1...n>] activities to [E;I] components respectively
        t; 
        dt = 1;
        tspan = [0 1];
        n = [2;2.2]; 
        k = 0.012;
        InSum;
        NullClines; % Handle - return nullclines for the specified input
        findIndex; % Handle for a function to find the values of this EI pair in the grouped column vector [<Exh 1;;;n>;<Inh 1;;;n>].
        % Cnorm = [0.5;0.5];
        % C50 = [50;50];
        % Cmax = [250;250];
        % Saturate = false;
    end
	
    methods 
        function objArr = Neurons(dimension, NStateVars, Wtdistfun, taus, ids) % Constructor (taus (seconds))
            if exist('dimension','var')
				objArr.Dim = dimension;
            end
            objArr.NNeuronGroups = prod(dimension); % EIpairArray must call the weight dist function here. Default Weight dist function can be pushed into EIpairArray class, along with tau and index functions
            if exist('NStateVars','var')
				objArr.NStateVarsPerGroup = NStateVars;
            end
            objArr.ids = [];
            if iscolumn(dimension)
                dimension = transpose(dimension);
            end
            if exist('ids','var')== 0
                objArr.ids = transpose(0:objArr.NNeuronGroups-1)*ones(1,length(dimension));
                cumdims = cumprod(dimension);
                tempid = (transpose(1:objArr.NNeuronGroups)*ones(1,size(dimension,2)) -1);
                tempid = floor(tempid / diag([1,cumdims(1:end-1)]));
                objArr.ids = 1 + round(tempid - floor(tempid/ diag(dimension))* diag(dimension));
            else
                objArr.ids = ids;
            end
            objArr.W = Wtdistfun();
%             size(objArr)
%             objArr.W = abs(objArr.W);
%             objArr.W(:,end/2+1 :end) = -objArr.W(:,end/2+1 :end);
%             size(objArr.W)
%             objArr.NNeuronGroups
            if exist('taus','var')==0
                taus = objArr.tau;
            end
            NNeuronGroups = objArr.NNeuronGroups;
            objArr.tau = zeros(NNeuronGroups*NStateVars,1);
            for i = 0:NStateVars-1
                objArr.tau(i*NNeuronGroups+1:(i+1)*NNeuronGroups) = taus(ceil((i*end+1)/NStateVars):ceil((i+1)*end/NStateVars));
            end          
        end
        
        
        function EI_ExtInp(objArr, Inp, tspan, R0, dt) % Inp = [2*NNeuronsx1] vector :: [Einputs(:); Iinputs(:)]
            objArr.Input = Inp;
            objArr.tspan = tspan;
            if exist('R0','var')==0
                R0 = objArr.R0;
            end
            objArr.R0 = zeros(objArr.NNeuronGroups*objArr.NStateVarsPerGroup,1);
            for i = 0:objArr.NStateVarsPerGroup-1
                objArr.R0(i*objArr.NNeuronGroups+1:(i+1)*objArr.NNeuronGroups) = R0(ceil((i*end+1)/objArr.NStateVarsPerGroup):ceil((i+1)*end/objArr.NStateVarsPerGroup));
            end
            if exist('dt','var')
                objArr.dt = dt; % in case of uniformly sampled inputs
            end
        end
        
        function SI = IpSumNonLinear(objArr ,y)   
            %y
            SI = max( objArr.W *y + objArr.Input , 0 ); 
            SI = objArr.k .*( SI .^ objArr.n );
%             if objArr.Saturate    
%                 SI = SI./(objArr.Cmax*0.01 + SI) .* objArr.Cmax;
%             end
        end
		
        function [ENC, INC] = NonLinearNC(obj, rectify, index )
            if exist('rectify','var') == 0
                rectify = true;
            end
            if exist('index','var')==0
                cumdim = cumprod(obj.Dim);
                index = ceil(obj.Dim*transpose([1,cumdim(1:end-1)])/2);
            end
            Wee = abs(obj.W(index,index));
            Wie = abs(obj.W(index+obj.NNeuronGroups,index));
            Wei = abs(obj.W(index,index+obj.NNeuronGroups));
            Wii = abs(obj.W(index+obj.NNeuronGroups,index+obj.NNeuronGroups));
            
            e = 2.^(-8:0.01:log(1000)/log(2));
            if rectify
                e = 2.^(-8:0.01:log(obj.Cmax(1))/log(2));
            end
            unse = e;
            if obj.Saturate 
                unse = e/obj.Cmax(1);
                unse = unse./(1-unse);
            end
            ENC = (-(unse/obj.k(1)).^(1/obj.n(1)) + Wee*e + obj.Input(1))/Wei;
            if rectify
                ENC = max(0,ENC);
            end
            ENC = [e;ENC];
            
            i = 2.^(-8:0.01:log(1000)/log(2));
            if rectify
                i = 2.^(-8:0.01:log(obj.Cmax(end))/log(2));
            end
            unsi = i;
            if obj.Saturate
                unsi = i/obj.Cmax(end);
                unsi = unsi./(1-unsi);
            end
            INC = ((unsi/obj.k(end)).^(1/obj.n(end)) + Wii*i - obj.Input(end))/Wie;
            if rectify
                INC = max(0,INC);
            end
            INC = [INC;i];
        end
        
%        function SI = IpSumNonLinearFullfield(objArr ,y, sigmapop)
%            %y
%            if isrow(objArr.Cnorm)
%                objArr.Cnorm = transpose(objArr.Cnorm);
%            end
%            if isrow(objArr.C50)
%                objArr.C50 = transpose(objArr.C50);
%            end
%            % Suppression - Full field stimulation - weights multiplied be
%            % a factor of 2pi/sigma^2
%            
%            sigma = ones(size(sigmapop))*(2*pi);
%            sigma = sigmapop.^2;
%%             objArr(:).W(1,:)*2*pi.*[sigma(1,1,1),sigma(1,1,end)]
%            SI = max( [objArr(:).W(1,:)*2*pi.*[sigma(1,1,1),sigma(1,1,end)];objArr(:).W(end,:)*2*pi.*[sigma(1,end,1),sigma(1,end,end)]] *y + objArr.Input , 0 ); 
%            SI = [objArr(:).k(1);objArr(:).k(end)] .* (SI.^[objArr(:).n(1);objArr(:).n(end)]);
%%             if objArr.Saturate
%%                 SI = SI./(1 + SI) .* objArr.Cmax;
%%             end
%        end
%		
%        function [ENC, INC] = NonLinearNCFullfield(obj, sigmapop, rectify)
%            
%            if exist('rectify','var') == 0
%                rectify = true;
%            end
%            if isrow(obj.Cnorm)
%                obj.Cnorm = transpose(obj.Cnorm);
%            end
%            if isrow(obj.C50)
%                obj.C50 = transpose(obj.C50);
%            end
%            sigma = ones(size(sigmapop))*(2*pi);
%            sigma = sigmapop.^2;
%            Wee = abs(obj.W(1,obj.ids*[1;transpose(obj.Dim(1:end-1))]));
%            Wei = abs(obj.W(1,obj.ids*[1;transpose(obj.Dim(1:end-1))] + obj.NNeuronGroups));
%            Wie = abs(obj.W(2,obj.ids*[1;transpose(obj.Dim(1:end-1))]));
%            Wii = abs(obj.W(2,obj.ids*[1;transpose(obj.Dim(1:end-1))] + obj.NNeuronGroups));
%            Wee = Wee*2*pi*sigma(1,1,1);
%            Wei = Wei*2*pi*sigma(1,1,2);
%            Wie = Wie*2*pi*sigma(1,2,1);
%            Wii = Wii*2*pi*sigma(1,2,2);
%            
%            e = 2.^(-8:0.01:log(1000)/log(2));
%            if rectify
%                e = 2.^(-8:0.01:log(obj.Cmax(1))/log(2));
%            end
%            unse = e;
%            if obj.Saturate
%                unse = e/obj.Cmax(1);
%                unse = unse./(1-unse);
%            end
%%             % Removing Normalization
%%             une= e/max(e);
%%             une = max((une./(1-une))*obj.C50(1)/obj.Cnorm(1),0);
%            ENC = ( -(unse/obj.k(1)).^(1/obj.n(1)) + Wee*e + obj.Input(1) )/Wei;
%            if rectify
%                ENC = max(0,ENC);
%            end
%            ENC = [e; ENC];
%            
%            i = 2.^(-8:0.01:log(1000)/log(2));
%            if rectify
%                i = 2.^(-8:0.01:log(obj.Cmax(end))/log(2));
%            end
%            unsi = i;
%            if obj.Saturate
%                unsi = i/obj.Cmax(end);
%                unsi = unsi./(1-unsi);
%            end
%%             % Removing Normalization
%%             uni = i/max(i);
%%             uni = max((uni./(1-uni)).*obj.C50(end)/obj.Cnorm(end), 0);
%            INC = ( (unsi/obj.k(end)).^(1/obj.n(end)) + Wii*i - obj.Input(end) )/Wie;
%            if rectify
%                INC = max(0,INC);
%            end
%            INC = [INC; i];
%            
%        end
%        
%        function SI = IpSumNonLinearPartial(objArr ,y, sigmapop, radius)
%            %y
%            if isrow(objArr.Cnorm)
%                objArr.Cnorm = transpose(objArr.Cnorm);
%            end
%            if isrow(objArr.C50)
%                objArr.C50 = transpose(objArr.C50);
%            end
%            % Suppression - Full field stimulation - weights multiplied be
%            % a factor of 2pi/sigma^2
%            
%            sigma = ones(size(sigmapop))*(2*pi);
%            sigma = sigmapop.^2 .* (1-exp(-0.5.*(radius./sigmapop).^2));
%%             objArr(:).W(1,:)*2*pi.*[sigma(1,1,1),sigma(1,1,end)]
%            SI = (max( [objArr(:).W(1,:)*2*pi.*[sigma(1,1,1),sigma(1,1,end)];objArr(:).W(end,:)*2*pi.*[sigma(1,end,1),sigma(1,end,end)]] *y +[objArr(:).Input(1); objArr(:).Input(end)] , zeros(size(y)) )); 
%            SI = [objArr(:).k(1);objArr(:).k(end)] .* (SI.^[objArr(:).n(1);objArr(:).n(end)]);
%            if objArr.Saturate
%                SI = SI./(1 + SI) .* objArr.Cmax;
%            end
%        end
%		
%        function [ENC, INC] = NonLinearNCPartial(obj, sigmapop,radius, rectify)
%            if exist('rectify','var') == 0
%                rectify = true;
%            end
%            if isrow(obj.Cnorm)
%                obj.Cnorm = transpose(obj.Cnorm);
%            end
%            if isrow(obj.C50)
%                obj.C50 = transpose(obj.C50);
%            end
%            sigma = ones(size(sigmapop))*(2*pi);
%            sigma = sigmapop.^2 .* (1-exp(-0.5.*(radius./sigmapop).^2));
%            
%            Wee = abs(obj.W(1,obj.ids*[1;transpose(obj.Dim(1:end-1))]));
%            Wei = abs(obj.W(1,obj.ids*[1;transpose(obj.Dim(1:end-1))] + obj.NNeuronGroups));
%            Wie = abs(obj.W(2,obj.ids*[1;transpose(obj.Dim(1:end-1))]));
%            Wii = abs(obj.W(2,obj.ids*[1;transpose(obj.Dim(1:end-1))] + obj.NNeuronGroups));
%            Wee = Wee*2*pi*sigma(1,1,1);
%            Wei = Wei*2*pi*sigma(1,1,2);
%            Wie = Wie*2*pi*sigma(1,2,1);
%            Wii = Wii*2*pi*sigma(1,2,2);
%            
%            e = 2.^(-8:0.01:log(1000)/log(2));
%            if rectify
%                e = 2.^(-8:0.01:log(obj.Cmax(1))/log(2));
%            end
%            unse = e;
%            if obj.Saturate
%                unse = e/obj.Cmax(1);
%                unse = unse./(1-unse);
%            end
%%             % Removing Normalization
%%             une= e/max(e);
%%             une = max((une./(1-une))*obj.C50(1)/obj.Cnorm(1),0);
%            ENC = ( -(unse/obj.k(1)).^(1/obj.n(1)) + Wee*e + obj.Input(1) )/Wei;
%            if rectify
%                ENC = max(0,ENC);
%            end
%            ENC = [e; ENC];
%            
%            i = 2.^(-8:0.01:log(1000)/log(2));
%            if rectify
%                i = 2.^(-8:0.01:log(obj.Cmax(end))/log(2));
%            end
%            unsi = i;
%            if obj.Saturate
%                unsi = i/obj.Cmax(end);
%                unsi = unsi./(1-unsi);
%            end
%%             % Removing Normalization
%%             uni = i/max(i);
%%             uni = max((uni./(1-uni)).*obj.C50(end)/obj.Cnorm(end), 0);
%            INC = ( (unsi/obj.k(end)).^(1/obj.n(end)) + Wii*i - obj.Input(end) )/Wie;
%            if rectify
%                INC = max(0,INC);
%            end
%            INC = [INC; i];
%            
%        end
%        
        function dR_dt = EI_update(objArr, t, y)  % 'y' is a [2*#EIpairs]x1 column vector of activities [<Exc1;2;;;>; <Inh1;2;;;>]
            y_latest = y(:,end);
%             size(y_latest),
%             y_latest(1), 
%             t,
            dR_dt = (-y_latest + ...
                objArr.InSum(objArr, ...
                 t(end), ...
                 y_latest )...
                 )./objArr.tau;
%             disp([dR_dt, Rlatest, [dR_dt(1).*objArr(:).tau(1);dR_dt(2).*objArr(:).tau(end)]+Rlatest])
        end
        
%        function ode(obj)
%            if isempty(obj.t) || obj.tspan(1)< obj.t(end)
%                obj.t = []; obj.R = [];
%            end
%            [ts, ys] = ode45(@(t,y) obj.EI_update(t,y), obj.tspan, obj.R0);
%            obj.t = [obj.t,transpose(ts)];
%            obj.R = [obj.R,transpose(ys)];
%        end
        
%        function collect(obj)  % Testing class references and debugging variables
%            if isempty(obj.R)
%                obj.E0 = obj.R0(1:length(R0)/2);
%                obj.I0 = obj.R0(length(R0)/2 + (1:length(R0)/2));
%                obj.E = obj.E0;
%                obj.I = obj.I0;
%                obj.t = [0];
%            end
%            obj.E = [obj.E,obj.E(:,end) + EI_update(obj,obj.t,obj.E).* obj.dt];
%            obj.I = [obj.I,obj.I(:,end) + EI_update(obj,obj.t,obj.I).* obj.dt];
%            obj.t = [obj.t,obj.t(end)+obj.dt];
%        end
        
%        function integrate(obj) % Testing External input summation functions
%            for i = 1:3
%                obj.collect();
%            end
%        end
    end
end