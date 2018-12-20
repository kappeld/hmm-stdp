function [net,Z,P] = snn_sample_pct( net, data, ct )
% Draws a spike from sem network - continuous time spiking hmm sems.
%
% [net,Z,P] = snn_sample_ct( net, data, ct )
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
%   tau           0.02         time constant of epsp window
%   lambda        100          network spike rate
%   refrac_time   0.000        absolute refractory time
%
% @fields:
%   At    zeros   [2,0]        log activity of network at spike times
%   

    t = data.time(ct(1));
    i = 1;
    
    time_range = data.time(ct(end))-data.time(ct(1));
    
    Delta = [];
    num_spikes = [];
    
    while isempty(num_spikes);
    
        Delta = [ Delta; exprnd( 1./(net.lambda), 100, 1 )+net.refrac_time ]; %#ok<AGROW>
    
        num_spikes = find( cumsum( Delta ) > time_range, 1, 'first' );
    end
    
    num_spikes = num_spikes-1;
    
    Delta = Delta(1:num_spikes);
        
    %num_spikes = ceil(time_range*net.lambda)-1;
    %Delta = 1/net.lambda*ones( 1, num_spikes );
    
    Z = zeros(2,num_spikes);
    P = zeros(net.num_neurons+1,num_spikes);
    net.At = zeros(2,num_spikes);
    net.Rt = zeros(net.num_samples+1,num_spikes);
    
    last_input_spike = t;
    
    hX = zeros( net.num_inputs, 1 );
    hX(1) = 1;
    
    hZ = zeros( net.num_neurons, 1 );
    
    eta = net.eta;
    
    W_nexp = exp( -net.W );
    V_nexp = exp( -net.V );
    V0_exp = exp( net.V0 );
    
    net.d_We = zeros(net.num_neurons,net.num_inputs);
    net.d_Wi = zeros(1,net.num_inputs);
    net.d_Ve =  zeros(net.num_neurons,net.num_neurons);
    net.d_Vi =  zeros(1,net.num_neurons);
    net.d_V0 = zeros(net.num_neurons,1);

    for j = 1:length(Delta)
        
        t = t+Delta(j);
        
        while (i < size(data.Xt,2)) && (t > data.Xt(2,i))
            hX = hX*exp(-(data.Xt(2,i)-last_input_spike)/net.tau);
            hX(data.Xt(1,i)) = hX(data.Xt(1,i)) + 1;
            last_input_spike = data.Xt(2,i);
            i = i+1;
        end
        
        hZ = hZ*exp(-Delta(j)/net.tau);
        hX = hX*exp(-(t-last_input_spike)/net.tau);
        
        last_input_spike = t;
        
        U = net.W*hX + net.V*hZ - net.V0;
        [P_t,A] = wta_softmax( U );
        
        Z_i = wta_draw(P_t);
        
        [k,v] = find(Z_i); %#ok<NASGU>
        
        Z(:,j) = [k;t];
        P(:,j) = [P_t;t];
        
        net.d_Wi = net.d_Wi + hX';
        net.d_Vi = net.d_Vi + hZ';
        
        net.d_We(k,:) = net.d_We(k,:) + min( W_nexp(k,:).*hX', 1/eta );
        net.d_Ve(k,:) = net.d_Ve(k,:) + min( V_nexp(k,:).*hZ', 1/eta );
        net.d_V0      = net.d_V0      + (Z_i-V0_exp)./max(eta,V0_exp);

        net.At(:,j) = [A,t];
        
        hZ(k) = hZ(k) + 1;
    end
end
