function [net,Z,P] = snn_sample_sempsp( net, data, ct )
% Draws a spike from sem network using psps.
%
% [net,Z,P] = snn_sample_sempsp( net, data, ct )
%
% Calculate membrane potentials net.U using post synaptic
% potentials with decay time constant net.tau. Each cell
% produces independently drawn spikes drawn from a Poisson
% process, which rate is proportional to the current
% output probability. The overall network spike rate is
% kept fixed to net.base_rate.
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
%   tau        0.02           membrane time constant
%   base_rate  100            network base fire rate
%   
%
% @fields:
%   U    zeros    [num_neurons,1]    membrane potentials
%   Y    zeros    [num_inputs,1]     synaptic activation
%   Psp  zeros    [num_neurons,1]    post synaptic potentials
%   Par  zeros    [num_neurons,1]    log partition function
%   A    zeros    [1,1]              network activation
%

    delta_t = 1/data.sample_rate;
    
    alpha = exp( -delta_t./net.tau );
    
    U = zeros(net.num_neurons,length(ct));
    
    for t = 1:length(ct)
        
        net.Psp = net.Psp*alpha;
        net.Y = net.Y*alpha;
        
        net.Psp = net.Psp + net.W*data.X(:,ct(t));
        net.Y = net.Y + data.X(:,ct(t));
        
        U(:,t) = net.Psp + net.W0 - net.Par;
    end
    
    [P,A] = wta_softmax(U);
    
    nextspike=exprnd(1./(net.base_rate*P+10));
    
    Z = nextspike < delta_t;
    
    net.A = A(:,end);
    net.U = U(:,end);
end
