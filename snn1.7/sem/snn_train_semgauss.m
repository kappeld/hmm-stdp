function [net,Z,P] = snn_train_semgauss( net, data, ct )
% Stochastic EM learning with variance tracking and bias weight update.
%
% [net,Z,P] = snn_train_semgauss( net, data, ct )
%
%
% @parameters:
%   eta                 0.05  learn rate
%   sample_method       'st'    default sample method
%   performance_method  ce      default performance method
%   temperature         1       initial weight temperature
%   train_w0            true    train bias weights
%
% @fields:
%   W         rand[temperature]  [num_neurons,num_inputs]  network weights
%   W0        rand[-1]           [num_neurons,1]           network bias
%   NumIter   zeros              [num_neurons,num_inputs]  # of iterations
%   S         zeros              [num_neurons,num_inputs]  weight means
%   Q         ones               [num_neurons,num_inputs]  weight variance
%   EtaW      zeros              [num_neurons,num_inputs]  weight learn-rates
%   S0        zeros              [num_neurons,1]           bias means
%   Q0        ones               [num_neurons,1]           bias variance
%   EtaW0     zeros              [num_neurons,1]           bias learn-rates
%   Par       zeros              [num_neurons,1]           log partition function
%

    num_samples = length(ct);

    if (nargin<3)
        ct = 1:size(data.Y,2);
    end

    P = zeros(net.num_neurons,num_samples);
    Z = zeros(net.num_neurons,num_samples);
    
    eta = net.EtaW;
    eta_W0 = net.EtaW0;
    
    for i=1:num_samples

        [net,Z_i,P(:,i)] = ...
            net.p_sample_fcn( net, data, ct(i) );
        
        Z(:,i) = Z_i;
        
        z_k = find(Z_i)';

        for k = z_k
            
            net.iteration = net.iteration+1;
            net.NumIter(k) = net.NumIter(k)+1;
            
            S_k = net.S(k,:);
            Q_k = net.Q(k,:);

            eta = net.eta*(Q_k-S_k.^2)./(exp(-S_k)+1);
            
            net.W(k,:) = net.W(k,:) + eta.*(net.Y'-net.W(k,:));

            net.S(k,:) = S_k + eta.*(net.W(k,:)-S_k);
            net.Q(k,:) = Q_k + eta.*(net.W(k,:).^2-Q_k);
        end

        net.Par = sum( net.W.^2, 2 )./2;
        
        if net.train_w0
            
            eta_w = net.eta*(net.Q0-net.S0.^2)./(exp(-net.S0)+1);

            eta_W0 = eta_w;

            net.W0 = net.W0 + eta_W0.*(Z_i-exp(net.W0))./max(eta_W0,exp(net.W0));

            net.S0 = net.S0 + eta_W0.*(net.W0-net.S0);
            net.Q0 = net.Q0 + eta_W0.*(net.W0.^2-net.Q0);
        end
    end
    
    net.EtaW = eta;
    net.EtaW0 = eta_W0;
end
