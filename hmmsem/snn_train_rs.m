function net = snn_train_rs( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems - rejection sampling.
%
% net = snn_train_ct( net, data, ct )
%
%
% @parameters:
%   eta                 0.005    learn rate
%   sample_method       'ct'     default sample method
%   performance_method  'ct'     default performance method
%   self_inhibition     0        fix value of self-recurrent weights
%   use_iw              true     use importance weight
%   resample            true     resample path if rejected
%   iw_track_speed      0.000001 speed of tracking importance weight mean
%   num_samples_min     10       min aspired number of samples
%   num_samples_max     10       max aspired number of samples
%   clip_weights        false    clip weights at -5
%   p_required_fields   [R,At,d_W,d_V,d_V0,SW_new,QW_new,SV_new,QV_new,S0_new,Q0_new]
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W        rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V        rand[-0.01]      [num_neurons,num_neurons] recurrent network weights
%   V0       rand[-0.01]    [num_neurons,1]           network bias
%   d_W      zeros          [num_neurons,num_inputs]  feedforward weight change
%   d_V      zeros          [num_neurons,num_neurons] recurrent weight change
%   d_V0     zeros          [num_neurons,1]           bias change
%   SW_new   zeros          [num_neurons,num_inputs]  new weight mean
%   QW_new   ones           [num_neurons,num_inputs]  new weight variance
%   SV_new   zeros          [num_neurons,num_neurons] new weight mean
%   QV_new   ones           [num_neurons,num_neurons] new weight variance
%   S0_new   zeros          [num_neurons,1]           new weight mean
%   Q0_new   ones           [num_neurons,1]           new weight variance
%   r_mean   rand[-3]       [1,1]                     current estimate of iw norm. (log)
%   biased_it       zeros   [0,0]                     iterations that were biased   
%   asp_num_samples rand[num_samples_min,num_samples_max] [1,1]  aspired number of resamples
%
%
% David Kappel
% 24.05.2011
%

    if ~isfield(net,'p_r_mean_init')
        net.p_r_mean_init = net.r_mean;
    end

    if net.use_iw && (net.asp_num_samples > 0)
        R_n = wta_draw( exp( double([ data(:).R ]') + net.r_mean ) );
        
        net.got_sample = any(R_n);

        if net.got_sample
            net.r_mean = net.r_mean - net.asp_num_samples*net.iw_track_speed;
            
            if ( sum(  exp( double([ data(:).R ]') + net.r_mean ) ) > 1 )
                warning('run is biased!');
                net.biased_it = [ net.biased_it, net.iteration+1 ];
            end

        else
            net.r_mean = net.r_mean + net.iw_track_speed;
        end
    else
        R_n = true(length(data),1);
        net.got_sample = true;
    end
    
    for i = find( R_n )'
        
        net.iteration = net.iteration + 1;

        if net.clip_weights
            net.W = max(-5,net.W + data(i).d_W);
            net.V = max(-5,net.V + data(i).d_V);
            net.V0 = max(-5,net.V0 + data(i).d_V0);
        else
            net.W = net.W + data(i).d_W;
            net.V = net.V + data(i).d_V;
            net.V0 = net.V0 + data(i).d_V0;
        end
                
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
    end    
end
