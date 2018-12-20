function [net,Z,P] = snn_sample_sp_inhibition_rec( net, data, ct )
% Spiking inhibition through recurrent spiknig neuron
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
%   I_rate            20             inhibitory base rate
%   refractory_period 0              refractory period
%   tau               0.02           membrane time constant
%   tau_I_r           0.0001         inhibitory rise time
%   tau_I_f           0.02           inhibitory fall time
%   w_fbI             3              feed back inh. weights
%   w_ffI             1              feed forward inh. weights
%   delay_I           1              inhibitory delay
%
% @fields:
%   U     zeros    [num_neurons,1]    Membrane potentials
%   A     zeros    [1,1]              Network activation
%   Y     zeros    [num_inputs,1]     Synaptic activation
%   I     zeros    [1,1]              Inh. membrane potential
%   ur_I  zeros    [1,1]              Inh. rinsing potential
%   uf_I  zeros    [1,1]              Inh. falling potential
%   z_I   zeros    [1,1]              Inh. output spikes
%   d_I   zeros    [1,delay_I]        Inhibitory delay line
%   t     zeros    [1,1]              Current simulation time step
%


    delta_t = 1/data.sample_rate;
    
    alpha = exp( -delta_t./net.tau );
    alpha_I_r = exp( -delta_t./net.tau_I_r );
    alpha_I_f = exp( -delta_t./net.tau_I_f );
    
    Z = zeros( net.num_neurons, length(ct) );
    P = zeros( net.num_neurons, length(ct) );
    
    for t = 1:length(ct)

        t_I = mod( net.t, net.delay_I )+1;
        
        % update z-neurons        
        net.Y = net.Y + data.X(:,ct(t));
        
        net.z_I = net.d_I(t_I);
        
        %if ( net.z_I > 0 )
        %    net.U(:) = 0;
        %end
        
        net.U = net.U + (0.5*net.W+0.5)*data.X(:,ct(t)) - net.z_I*net.w_fbI;
        
        %nextspike = exprnd(1./(net.base_rate*exp(net.U)));
        nextspike = exprnd(1./(net.base_rate*net.U));

        Z(:,t) = nextspike < delta_t;
        P(:,t) = exp(net.U);
        
        % update inhibitory neuron
        net.I = net.uf_I - net.ur_I;
        
        nextspike = exprnd(1./(net.I_rate*exp(net.I)));
        
        net.d_I(t_I) = (nextspike < delta_t);
        
        net.ur_I = net.ur_I + sum(Z(:,t))*net.w_ffI;
        net.uf_I = net.uf_I + sum(Z(:,t))*net.w_ffI;
        
        net.ur_I = net.ur_I*alpha_I_r;
        net.uf_I = net.uf_I*alpha_I_f;
        net.U = net.U*alpha;
        net.Y = net.Y*alpha;
        
        net.t = net.t + 1;
    end
    
    net.A = sum(P(:,end));

end
