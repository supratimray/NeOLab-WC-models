%% Generate ISN model
function pop = ISN_WCRinzel2017(nsimulations, argWeights, gamma, k, ks, theta, thetas)
    if ~exist('nsimulations','var')
        nsimulations = 1%[1, 1];
    end
    if exist('argWeights','var')==0
        Weights = zeros([4,4]);
        Weights(1:2,3:4) = [4, -6 ; 3, -3]; % 2-variable form - Keeley NYU Dissertation
        %[3.5, -5 ; 3.5, -3];%Rinzel 2019 Firing Rate Models paper.
        %[3, -3 ; 1, -0]%Rinzel 2017 dual gamma model paper. 
        %[16, -12 ; 15, -3];%WC1972 
        % [16, -26 ; 20, -1]; %JS2014
    else
        if all(size(argWeights) == [2,2])
            Weights = zeros([4,4]);
            Weights(1:2,3:4) = argWeights;
        else
            Weights = argWeights;
        end
        Weights
    end
    pop = GenerateNeuronPopulation; % EI_populationArray;
    % TODO: configure ssnpop for experiment
    pop.connectivityMaxW = Weights;
    pop.connectivitySigma = 0; % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    taus = [0, 0, 5, 10]*1e-3;
    nStateVars = 4; % corresponding to E and I cells' rate variables, followed by their respective Synaptic variables: r_E, r_I, s_E, s_I
    pop.setupEIPopulation(nsimulations, nStateVars, 'Keeley 2017', taus);
    function op = RinzelSynapticSum(objArr, x, gamma, theta, thetas, k, ks)%, s0)
        %     x
        
        OPsize = size(x);
        xsum = max( objArr.W *x + objArr.Input , 0 );
        x = reshape(x, [length(x)/4, 4]);
        xsum = reshape(xsum, [length(xsum)/4, 4]);
        xsum = xsum(:,1:2);
        rop = 1./(1+exp((theta-xsum)*diag(1./k)));

%         r = x(:,1:2); 
        s = x(:,3:4);
        gamma = reshape(gamma, [1, 2]);
        % s0= reshape(s0, [1, 2]);
        sop = (1-s)./(1+exp((thetas-rop)*diag(1./ks)))*diag(gamma);
        % sop = sop +s0;
        op = reshape([rop,sop], OPsize);
    end

%     function op = TestJSRinzelSynapticSum(objArr, x, gamma, theta, thetas, k, ks)%, s0)
%         %     x
%         
%         OPsize = size(x);
%         xsum = max( objArr.W *x + objArr.Input , 0 );
%         x = reshape(x, [length(x)/4, 4]);
%         xsum = reshape(xsum, [length(xsum)/4, 4]);
%         xsum = xsum(:,1:2);
%         rop = xsum;%1./(1+exp((theta-xsum)));
% 
%         r = x(:,1:2); s = x(:,3:4);
%         gamma = reshape(gamma, [1, 2]);
%         % s0= reshape(s0, [1, 2]);
%         sop = 1./(1+exp([5,20]-r)); %(1-s)./(1+exp((thetas-r)*diag(1./ks)))*diag(gamma);
%         % sop = sop +s0;
%         op = reshape([rop,sop], OPsize);
%     end

    if ~exist('gamma','var')
        gamma = [3,4]; 
    end
    if ~exist('k','var')
        k = [0.1 0.1]; 
    end
    if ~exist('ks','var')
        ks=[0.1 0.1];
    end
    if ~exist('theta','var')
        theta = [0.2 0.4]; 
    end
    if ~exist('thetas','var')
        thetas = [0.9, 0.9];%s0 = [0,0];
    end
    %pop.Population_Name = [pop.Population_Name, ' _Gamma=',num2str(gamma),' _k=', num2str([k,ks]),' _theta=', num2str([theta,thetas])];%' _s0=', num2str(s0)];
    pop.EIpairs.InSum = @( objArr, t, y) RinzelSynapticSum( objArr, y, gamma, theta, thetas, k, ks)%, s0);
end