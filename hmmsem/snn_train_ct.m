function net = snn_train_ct( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems.
%
% net = snn_train_ct( net, data, ct )
%
%
% @parameters:
%   eta                 0.005   learn rate
%   sample_method       'ct'    default sample method
%   performance_method  'ct'    default performance method
%   self_inhibition     0       fix value of self-recurrent weights
%   use_iw              true    use importance weights
%   p_required_fields   [R,d_W,d_V,d_V0]
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W     rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V     zeros          [num_neurons,num_neurons] recurrent network weights
%   V0    rand[-0.01]    [num_neurons,1]           network bias
%   d_W   zeros          [num_neurons,num_inputs]  feedforward weight change
%   d_V   zeros          [num_neurons,num_neurons] recurrent weight change
%   d_V0  zeros          [num_neurons,1]           bias change
%

    if net.use_iw
        R_all = [ data(:).R ];
        R_all = mk_normalised( exp(R_all-max(R_all)) );
    else
        R_all = mk_normalised( ones(length(data),1) );
    end
    
    net.iteration = net.iteration + 1;
    
    eta = min( net.eta, (1/net.iteration) );

    for i=ct
        net.W = net.W + eta*R_all(i)*data(i).d_W;
        net.V = net.V + eta*R_all(i)*data(i).d_V;
        net.V0 = net.V0 + eta*R_all(i)*data(i).d_V0;
    end
    
    if ~isnan( net.self_inhibition )
        net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
    end
    
    net.got_sample = true;
    
end
