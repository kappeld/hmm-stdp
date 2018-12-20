function [net,Z,P] = snn_sample_sempoismtl( net, data, ct )
% Draws a spike assuming poission distributed inputs with multiple time lags.
%
% [net,Z,P] = snn_sample_sempoismtl( net, data, ct )
%
% Generates the output of a given wta-network by drawing
% a winner neuron from the current output probabilities
% assuming poisson distributed inputs.
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
%   lead_in        20             synaptic integration length
%   base_rate      100            network base fire rate
%   num_time_lags  3              number of time lags
%
% @fields:
%   U   zeros   [num_neurons,1]              Membrane potentials
%   A   zeros   [1,1]                        Network activation
%   Y   zeros   [num_inputs,1]               synaptic activation
%

    delta_t = 1/data.sample_rate;
    
    f_size = net.lead_in/net.num_time_lags;
    
    num_syn = net.num_inputs/net.num_time_lags;
        
    U = zeros( net.num_neurons, length(ct) );
    
    for t = 1:length(ct)

        for l = 0:(net.num_time_lags-1)
            net.Y(l*num_syn+(1:num_syn)) = sum( data.X(:,ct(t)-(0:(f_size-1))-(f_size*l)), 2 );
        end
        U(:,t) = net.W*net.Y + net.W0 - sum( exp(net.W), 2 );
    end
    
    [P,A] = wta_softmax(U);
    
    nextspike=exprnd(1./(net.base_rate*P));
    
    Z = nextspike < delta_t;
    
    net.U = U(:,end);
    net.A = A(:,end);

end
