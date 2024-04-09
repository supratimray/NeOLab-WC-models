%% Generate ISN model
function pop = ISN_CRS(nsimulations, LocalWeights, LateralWeights, InputWeights)
    if ~exist('nsimulations','var')
        nsimulations = [2,2];
    end
    if exist('LocalWeights','var')==0
        LocalWeights = [20, -30, -20; 15, -12, -10; 0, -0, -0];
    end
    if exist('LateralWeights','var')==0
        LateralWeights = [50; 100; 0];
		LateralWeights = LateralWeights*[1,-1,0];
    end
	Weights = LocalWeights + LateralWeights;
    pop = GenerateNeuronPopulation; % EI_populationArray;
    % TODO: configure ssnpop for experiment
    pop.connectivityMaxW = Weights;
    pop.connectivitySigma = 0; % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    taus = [8,12,12]*1e-3;
    nStateVars = 3; % corresponding to E, I and S cells' state variables
    pop.setupEIPopulation(nsimulations, nStateVars, 'VinayCRS model',  taus);

    
    if exist('InputWeights','var')==0
        InputWeights = eye(3);%specifying IE, II and IG directly;  
		% InputWeights = [Weff, Wefb, 0; Wiff, Wifb, 0; 0, 0, Wsfb]% full model
	end
	
	function op = sigmoid(objArr, y, m, theta, InputWeights)
		% y
		OPsize = size(y);
		input = reshape(objArr.Input, [objArr.NNeuronGroups, objArr.NStateVarsPerGroup])*transpose(InputWeights);
		x = max( objArr.W*y + input(:), 0);
        
        x = reshape(x, [objArr.NNeuronGroups, objArr.NStateVarsPerGroup]);
        theta = reshape(theta, [1, objArr.NStateVarsPerGroup]);
        op = x-theta;
        m = reshape(m, [1, objArr.NStateVarsPerGroup]);
        op = op*diag(m);
        op = 1./(1+exp(-op));
        op = reshape(op, OPsize);
    end
    m = [1, 1.2, 1.5]; theta = [5, 20, 25];
    
    % pop.Population_Name = [pop.Population_Name, '_m=', num2str(m),' _Theta=',num2str(theta)];
    pop.EIpairs.InSum = @( objArr, t, y) sigmoid(objArr, y, m, theta, InputWeights);
end