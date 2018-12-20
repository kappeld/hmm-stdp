function net = snn_train_ct_tiw( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems, iw tracking.
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
%   iw_track_speed      0.005   speed of tracking importance weights (0..1)
%   p_required_fields   [R,d_W,d_V,d_V0]
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W      rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V      rand[-0.01]    [num_neurons,num_neurons] recurrent network weights
%   V0     rand[-0.01]    [num_neurons,1]           network bias
%   d_W    zeros          [num_neurons,num_inputs]  feedforward weight change
%   d_V    zeros          [num_neurons,num_neurons] recurrent weight change
%   d_V0   zeros          [num_neurons,1]           bias change
%   r_mean zeros          [1,0]                     current estimate of iw normalisation
%

    if net.use_iw
        R_all = [ data(:).R ];
        R_n =  zeros( length(data),1 );
        for i=ct
           if isempty(net.r_mean)
               net.r_mean = mean( R_all(i) );
           end
           net.r_mean = (1-net.iw_track_speed)*net.r_mean + net.iw_track_speed*R_all(i);
           R_n(i) = exp( R_all(i) - net.r_mean )/length(data);
        end
%         for i=ct
%             R_max = max( R_all );
%             
%             if isempty(net.r_mean)
%                 net.r_mean = mean( R_all(i) );
%             end
%             net.r_mean = log( (1-net.iw_track_speed)*exp( net.r_mean - R_max ) + ...
%                               net.iw_track_speed*exp( R_all(i) - R_max ) ) + R_max;
%             R_n(i) = exp( R_all(i) - net.r_mean );
%         end
        
        R_n = min(1,R_n);
    else
        R_n = mk_normalised( ones(length(data),1) );
    end
    
    net.iteration = net.iteration + 1;
    
    eta = min( net.eta, (1/net.iteration) );

    for i=ct
        net.W = net.W + eta*R_n(i)*data(i).d_W;
        net.V = net.V + eta*R_n(i)*data(i).d_V;
        net.V0 = net.V0 + eta*R_n(i)*data(i).d_V0;
    end
    
    if ~isnan( net.self_inhibition )
        net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
    end
    
    net.IW = R_n;
end
