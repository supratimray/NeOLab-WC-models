%% Generate ISN model
function pop = ISN_JXK(nsimulations, Weights, MN, r)
    if ~exist('nsimulations','var')
        nsimulations = [5,5];
    end
    Wflag = false;
    if exist('Weights','var')==0
        Wflag = true;
    elseif isempty(Weights)
        Wflag = true;
    end
    if Wflag
        Weights = [1.5, -3.25, 0.25; 3.5, -2.5, 0.5; 0.6, -0, 0];
    end
    pop = GenerateNeuronPopulation; % EI_populationArray;
    % TODO: configure ssnpop for experiment
    pop.connectivityMaxW = Weights;
    pop.connectivitySigma = 0; % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    taus = [6,15,19]*1e-3;
    nStateVars = 3; % corresponding to E, I and S cells' state variables
    pop.setupEIPopulation(nsimulations, nStateVars, 'Jia Xing Kohn 2013',  taus);

    function op = JXK2013summation(objArr, t, x, MN, r)
        OPsize = size(x);
        Ne = 1; Ni = 2; Ng = 3; % State variable indices corresponding to E and I cell rates
        NNeuronGroups = objArr.NNeuronGroups;
        NStateVarsPerGroup = objArr.NStateVarsPerGroup;

        x = max( x, 0); % H(x) in eqn(4) of the paper
        x = reshape(x, [NNeuronGroups, NStateVarsPerGroup]);
        x(:,1) = x(:,1).*(1-MN);
        poissoninputs = poissrnd(objArr.Input);

        x = objArr.W*x(:) + poissoninputs;  
        x = reshape(x, [NNeuronGroups, NStateVarsPerGroup]);
        x(:,3) = x(:,3).*(r.^2);
        op = reshape( x, OPsize);
    end


    if ~exist('MN','var')
        MN = 0; 
    end
    if ~exist('r','var')
        r = 5;
    end
    
    % pop.Population_Name = [pop.Population_Name, '_MN=', num2str(MN),'_r=',num2str(r)];
    pop.EIpairs.InSum = @( objArr, t, y) JXK2013summation(objArr, t, y, MN, r);
end