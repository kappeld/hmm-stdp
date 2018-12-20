function net = snn_train_rsw0( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems - rejection sampling.
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
%   iw_track_speed      0.05   speed of tracking importance weight mean
%   p_required_fields   [R,d_W,d_V,d_V0,d_W0,SW_new,QW_new,SV_new,QV_new,S0_new,Q0_new]
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W        rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V        rand[-0.01]    [num_neurons,num_neurons] recurrent network weights
%   V0       rand[-0.01]    [num_neurons,1]           network bias
%   W0       rand[-0.01]    [num_neurons,1]           network bias
%   d_W      zeros          [num_neurons,num_inputs]  feedforward weight change
%   d_V      zeros          [num_neurons,num_neurons] recurrent weight change
%   d_V0     zeros          [num_neurons,1]           bias change
%   SW_new   zeros          [num_neurons,num_inputs]  new weight mean
%   QW_new   ones           [num_neurons,num_inputs]  new weight variance
%   SV_new   zeros          [num_neurons,num_neurons] new weight mean
%   QV_new   ones           [num_neurons,num_neurons] new weight variance
%   S0_new   zeros          [num_neurons,1]           new weight mean
%   Q0_new   ones           [num_neurons,1]           new weight variance
%   r_mean   zeros          [1,1]                     current estimate of iw normalisation
%
%
% David Kappel
% 24.05.2011
%
%

    iw_track_speed = net.iw_track_speed;
    
    %iw_track_speed = 0.1/(net.iteration + 1);

    if net.use_iw
        
        R_all = [ data(:).R ];
        R_n =  zeros( length(data),1 );
        
        for i=ct

            p_accept = exp( double(R_all(i)) - net.r_mean );
            
            %iw_track_speed*exp(R_all(i));
            
            net.r_mean = (1-iw_track_speed)*net.r_mean + iw_track_speed*R_all(i);

            %net.r_mean = net.r_mean + log1p( iw_track_speed*expm1( double(R_all(i)) - net.r_mean ) );
        
            if ( p_accept > rand() )
                R_n(i) = 1;
                break;
            end
        end
        
    else
        R_n = ones(length(data),1);
    end
    
    net.got_sample = false;
    
    for i = find( R_n > 0 )'
        
        net.iteration = net.iteration + 1;

        net.W = net.W + data(i).d_W;
        net.V = net.V + data(i).d_V;
        net.V0 = net.V0 + data(i).d_V0;
        net.W0 = net.W0 + data(i).d_W0;
        
        if ~isnan( net.self_inhibition )
            net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
        end
        
        if net.use_variance_tracking
            net.SW = data(i).SW_new;
            net.QW = data(i).QW_new;
            net.SV = data(i).SV_new;
            net.QV = data(i).QV_new;
            net.S0 = data(i).S0_new;
            net.Q0 = data(i).Q0_new;
        end

        net.IW = R_n;
        net.got_sample = true;
    end
end
