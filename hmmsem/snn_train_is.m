function net = snn_train_is( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems - importance sampling.
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
%   iw_track_speed      0.05    speed of tracking importance weight mean
%   num_samples         10      number of samples for importance sampling
%   p_required_fields   [R,d_W,d_V,d_V0,SW_new,QW_new,SV_new,QV_new,S0_new,Q0_new]
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W        rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V        rand[-0.01]    [num_neurons,num_neurons] recurrent network weights
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
%   r_mean   zeros          [1,1]                     current estimate of iw normalisation
%
% David Kappel
% 24.05.2011
%

    R_all = wta_softmax( [ data(:).R ]' );
    
    net.r_mean = mean( [ data(:).R ], 2 );
    
    net.iteration = net.iteration + 1;

    net.IW = R_all;
    net.got_sample = true;
    
    if ~net.use_iw
        R_all = wta_draw( R_all );
    end
    
    for i = 1:length( R_all )
        net.W = net.W + R_all(i)*data(i).d_W;
        net.V = net.V + R_all(i)*data(i).d_V;
        net.V0 = net.V0 + R_all(i)*data(i).d_V0;
    end
    
    if net.use_variance_tracking
        
        SW = zeros(net.num_neurons,net.num_inputs);
        QW = zeros(net.num_neurons,net.num_inputs);
        SV = zeros(net.num_neurons,net.num_neurons);
        QV = zeros(net.num_neurons,net.num_neurons);
        S0 = zeros(net.num_neurons,1);
        Q0 = zeros(net.num_neurons,1);

        for i = 1:length( R_all )
            SW = SW + R_all(i)*data(i).SW_new;
            QW = QW + R_all(i)*data(i).QW_new;
            SV = SV + R_all(i)*data(i).SV_new;
            QV = QV + R_all(i)*data(i).QV_new;
            S0 = S0 + R_all(i)*data(i).S0_new;
            Q0 = Q0 + R_all(i)*data(i).Q0_new;
        end

        net.SW = SW;
        net.QW = QW;
        net.SV = SV;
        net.QV = QV;
        net.S0 = S0;
        net.Q0 = Q0;
    end
    
    if ~isnan( net.self_inhibition )
        net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
    end
end
