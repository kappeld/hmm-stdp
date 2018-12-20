function [net,Z,P] = snn_train_pct( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems.
%
% [net,Z,P] = snn_train_ct( net, data, ct )
%
%
% @parameters:
%   eta                 0.005   learn rate
%   sample_method       'ct'    default sample method
%   performance_method  'ce'    default performance method
%   num_samples         100     number of sampled paths
%   p_requirements      [batch_mode,continuous_time]
%
% @fields:
%   W     rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V     rand[-0.01]    [num_neurons,num_neurons] recurrent network weights
%   V0    rand[-0.01]    [num_neurons,1]           network bias
%   d_We  zeros   [num_neurons,num_inputs]   exitatory feedforward weight update
%   d_Wi  zeros   [1,num_inputs]             inhibitory feedforward weight update
%   d_Ve  zeros   [num_neurons,num_neurons]  exitatory recurrent weight update
%   d_Vi  zeros   [1,num_neurons]            inhibitory recurrent weight update
%   d_V0  zeros   [num_neurons,1]            bias update
%   R     zeros   [num_samples,1]            importance weights
%

    R = zeros(1, net.num_samples);

    d_W =  zeros(net.num_neurons,net.num_inputs,net.num_samples);
    d_V =  zeros(net.num_neurons,net.num_neurons,net.num_samples);
    d_V0 = zeros(net.num_neurons,1,net.num_samples);    
    net.R = zeros(net.num_samples,size(data.time,2));
    
    for i=1:net.num_samples
        
        [net,Z,P] = ...
            net.p_sample_fcn( net, data, ct );
        
        d_W(:,:,i) = net.d_We - repmat(net.d_Wi/size(P,2),net.num_neurons,1);
        d_V(:,:,i) = net.d_Ve - repmat(net.d_Vi/size(P,2),net.num_neurons,1);
        d_V0(:,i) = net.d_V0;
        
        R_all = cumsum( net.At(1,:) );
        
        R(i) = R_all(end);
        
        for j=1:size(net.At,2)
            net.R(i,ceil(net.At(2,j)*data.sample_rate):end) = R_all(j);
        end
    end

    R = mk_normalised( exp(R-max(R)) );
    
    net.iteration = net.iteration + 1;
    
    eta = net.eta;

    for i=1:net.num_samples
        net.W = net.W + eta*R(i)*d_W(:,:,i);
        net.V = net.V + eta*R(i)*d_V(:,:,i);
        net.V0 = net.V0 + eta*R(i)*d_V0(:,i);
    end
    
    %net.W = log( mk_stochastic( exp(net.W)' ) )';
    %net.V = log( mk_stochastic( exp(net.V)' ) )';
    %net.V0 = log( mk_normalised( exp(net.V0) ) );

    if ( net.self_inhibition ~= 0 )
        net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
    end
    
end
