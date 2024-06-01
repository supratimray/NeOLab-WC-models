%% Initiate and Configure JS model
function pop = ISN_JS2014(nsimulations, noisefile, taus, theta,Weights, A)
    if ~exist('nsimulations','var')
        nsimulations = 1;
    end
    if exist('Weights','var')==0
        Weights = [16, -26 ; 20, -1]; 
    end

    if ~exist('A','var')
        A = [1,1]; 
    end
    if ~exist('theta','var')
        theta = [5, 20];
    end

    if ~exist('noisefile','var')
        noisefile = [];
    end
    

    if isempty(noisefile)
        getnoise = @(t) 0;
    else
        % Load input/noise timeseries from file
        try
            load(noisefile,'noise','noise_t');
        catch
            try

                load(noisefile,'noise','noise_t');
            catch
                try

                    load(noisefile,'noise','noise_t');
                catch

                    load(noisefile,'noise','noise_t');
                end
            end
        end
        getnoise = @(t) noise(:,noise_t==t);
    end

    pop = GenerateNeuronPopulation; 
    pop.Population_Name = 'JS2014';
    pop.connectivityMaxW = Weights;
    pop.connectivitySigma = 0; % No interconnection between adjacent EI_pairs - multiple independent JS simulations
    % Makes each EI_pair connect only within itself (Each EI_pair can be a separate simulation).
    
    if ~exist('taus','var')
        taus = ([20,10]*1e-3); % seconds
    end
    
    % Setup neuron class object to record simulation outputs and states
    nStateVars = 2; % Firing rates of E and I population
    pop.setupEIPopulation(nsimulations, nStateVars, 'JS2014', taus);
    
    % Defining activation function to be used to update states of the
    % dynamical system
    function op = normadd_sigmoid(x, A, theta)
        x = x(:); 
        A = 1./A(:)'; %A = reshape(1./A, [1, 2]);
        op = [x(1:end/2), x(end/2+1:end)]-ones(numel(x)/2,1)*theta(:)';
        op2 = (-ones(numel(x)/2,1)*theta(:)');        

        op = op*diag(A);
	    op2 = op2*diag(A);
        
        op = 1./(1+exp(-op));
        op2 = 1./(1+exp(-op2));
        
        op = op(:) - op2(:); 
    end

    % Specify the Nonlinear summation part of the Activation function 'f' :
        % where dy/dt = tau*( -y + f(y,ip,t) ) or y = f(y,ip,t) ... latter
        % case is implemented for variables with tau=inf

    pop.EIpairs.InSum = @( objArr, t, y, recordflag) ...
        normadd_sigmoid( objArr.W *y ...
        + (objArr.Input) + getnoise(t)...
        , A, theta);
end