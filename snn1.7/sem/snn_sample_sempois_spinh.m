function [net,Z,P] = snn_sample_sempois_spinh( net, data, ct )
% Draws a spike from sem network assuming poission distributed inputs with refractory time.
%
% [net,Z,P] = snn_sample_sempois( net, data, ct )
%
% Generates the output of a given wta-network by drawing
% a winner neuron from the current output probabilities
% (standard version).
%
% inputs:
%   net:  A wta-network, see: wta-new()
%   data: A data structure to be simulated.
%   ct:   Current time index.
%
% output:
%   net:      The (modified) network structure
%   Z:        Network output spikes
%
% @parameters:
%   lead_in    30             synaptic integration length
%   base_rate  1000           network base fire rate
%   T_refrac   30             frefractory time (in time steps)
%
% @fields:
%   U    zeros    [num_neurons,1]    Membrane potentials
%   A    zeros    [1,1]              Network activation
%   Y    zeros    [num_inputs,1]     synaptic activation
%   Ts   zeros    [num_neurons,1]    last spike time
%   ct   zeros    [1,1]              current time
%

    delta_t = 1/data.sample_rate;
    
    f_size = net.lead_in;
    
    U = zeros( net.num_neurons, length(ct) );
    P = zeros( net.num_neurons, length(ct) );
    Z = false( net.num_neurons, length(ct) );
    
    for t = 1:length(ct)

        net.Y = sum( data.X(:,ct(t)-(0:(f_size-1))), 2 );
        U(:,t) = net.W*net.Y + net.W0; % - 0.5*sum( exp(net.W), 2 );
        
        net.ct = net.ct + 1;
        
        %U_exp = exp(0.1*U(:,t));
        %P(:,t) = U_exp./(1+U_exp);
        
        P(:,t) = exp(U(:,t));   
    
        nextspike=exprnd(1./(net.base_rate*P(:,t)));
    
        Z(:,t) = nextspike < delta_t;
    
        for l = find(Z(:,t))'
            net.Ts(l) = net.ct;
        end
    end
    
    net.U = U(:,end);
    net.A = sum( P(:,end) );

end
