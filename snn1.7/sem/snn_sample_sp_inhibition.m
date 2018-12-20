function [net,Z,P] = snn_sample_sp_inhibition( net, data, ct )
% Spiking lateral inhibition SEM sampling.
%
% [net,Z,P] = snn_sample_sp_inhibition( net, data, ct )
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
%   base_rate         20             network base fire rate
%   refractory_period 0              refractory period
%   tau               0.02           membrane time constant
%   alpha_I           0.98           inhibitory decay
%   w_I               3              inhibitory weight
%   delay_I           1              inhibitory delay
%
% @fields:
%   U    zeros    [num_neurons,1]    Membrane potentials
%   V    zeros    [num_neurons,1]    Normalised membrane potentials
%   Y    zeros    [num_inputs,1]     Synaptic activation
%   I    zeros    [1,1]              Inhibitory current
%   A    zeros    [1,1]              Network activation
%   d_I  zeros    [1,delay_I]        Inhibitory delay line
%


    delta_t = 1/data.sample_rate;
    
    alpha = exp( -delta_t./net.tau );
    
    Z = zeros( net.num_neurons, length(ct) );
    P = zeros( net.num_neurons, length(ct) );
    
    for t = 1:length(ct)
        
        net.Y = net.Y + data.X(:,ct(t));
        
        net.U = (0.5*net.W+0.5)*net.Y;
        %net.U = net.W*net.Y + net.W0;
        
        net.V = net.U - net.I;
        
        net.I = net.I*net.alpha_I;
        
        nextspike = exprnd(1./(net.base_rate*exp(net.V)));
        %nextspike = exprnd(1./(net.base_rate*net.V));

        Z(:,t) = nextspike < (1/data.sample_rate);
        P(:,t) = exp(net.V);

        t_I = mod( ct(t), net.delay_I )+1;

        net.I =  net.I + net.d_I(t_I);        
        net.d_I(t_I) = sum(Z(:,t))*net.w_I;
        
        net.Y = net.Y*alpha;
    end
    
    net.A = sum(P(:,end));

end
