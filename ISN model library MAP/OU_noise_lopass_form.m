%% Generate ISN model
function pop = OU_noise_lopass_form(nsimulations, noisefile, theta)
    if ~exist('nsimulations','var')
        nsimulations = 1; %[1, 1];
    end
    
    if ~exist('theta','var')
        theta = [1 1];
    end
    
    if ~exist('noisefile','var')
        getmeanip = @(t) 0;
        getnoise = @(t) 0;
    else
        noise = noisefile{1};
        meanip = noisefile{2};
        noise_t = noisefile{3};
        
        getmeanip = @(t) meanip(:,noise_t==t);
        getnoise = @(t) noise(:,noise_t==t);
    end

    nStateVars = 2; % corresponding to E and I cells' state variables
    
    pop = GenerateNeuronPopulation; % EI_populationArray;
    % TODO: configure ssnpop for experiment
    pop.connectivityMaxW = zeros(nStateVars);
    pop.connectivitySigma = 0; % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    if ~exist('taus','var')
        taus = 1./theta; 
    end

    pop.setupEIPopulation(nsimulations, nStateVars, 'OUpop', taus);

    function op = ou_filtForm(meanip, ip)
        op = meanip + ip;
        op = op(:);
    end
    pop.EIpairs.InSum = @( objArr, t, y, recordflag) ou_filtForm(getmeanip(t),getnoise(t));
end