%% Generate ISN model
function pop = ISN_WCRinzel2019(nsimulations, Weights, gamma, s0)
    if ~exist('nsimulations','var')
        nsimulations = 1%[1, 1];
    end
    if exist('Weights','var')==0
        Weights = zeros([4,4]);
        Weights(3:4,1:2) = [10, -3 ; 1, -3];%Rinzel 2017 dual gamma model paper. %[16, -12 ; 15, -3];%WC1972 % [16, -26 ; 20, -1]; %JS2014
    end
    pop = GenerateNeuronPopulation; % EI_populationArray;
    % TODO: configure ssnpop for experiment
    pop.connectivityMaxW = Weights;
    pop.connectivitySigma = 0; % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    taus = [5, 5, 3, 10]*1e-3;
    nStateVars = 4; % corresponding to E and I cells' rate variables, followed by their respective Synaptic variables: r_E, r_I, s_E, s_I
    pop.setupEIPopulation(nsimulations, nStateVars, 'Keeley 2019', taus);
    function op = sigmoid(objArr, x, gamma, s0)
        %     x
        OPsize = size(x);
        xsum = max( objArr.W *x + objArr.Input , 0 );
        x = reshape(x, [length(x)/4, 4]);
        xsum = reshape(xsum, [length(xsum)/4, 4]);
        xsum = xsum(:,3:4);
        r = x(:,1:2); s = x(:,3:4);
        gamma = reshape(gamma, [1, 2]);
        s0= reshape(s0, [1, 2]);
        sop = (r.*(1-s))*diag(gamma)+s0;
        rop = 1./(1+exp(-xsum));
        op = reshape([rop,sop], OPsize);
    end
    if ~exist('gamma','var')
        gamma = [1,1]; 
    end
    if ~exist('s0','var')
        s0 = [0,0];
    end
    %pop.Population_Name = [pop.Population_Name, ' _Gamma=',num2str(gamma),' _s0=', num2str(s0)];
    pop.EIpairs.InSum = @( objArr, t, y) sigmoid( objArr, y, gamma, s0);
end