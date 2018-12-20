function [net,Z,P] = snn_train_stw0( net, data, ct )
% Stochastic EM learning with variance tracking and bias weight update.
%
% [net,Z,P] = snn_update_st( net, data, ct )
%
%
% @parameters:
%   eta                 0.05    learn rate
%   sample_method       'st'    default sample method
%   performance_method  ce      default performance method
%   temperature         -0.01   initial weight temperature
%
% @fields:
%   W         rand[temperature]  [num_neurons,num_inputs]  network weights
%   W0        rand[temperature]  [num_neurons,1]           network bias
%   num_iter  zeros              [num_neurons,num_inputs]  # of iterations
%   S         zeros              [num_neurons,num_inputs]  weight means
%   Q         ones               [num_neurons,num_inputs]  weight variance
%   etaW      zeros              [num_neurons,num_inputs]  weight learn-rates
%   S0        zeros              [num_neurons,1]           bias means
%   Q0        ones               [num_neurons,1]           bias variance
%   etaW0     zeros              [num_neurons,1]           bias learn-rates
%

    num_samples = length(ct);

    if (nargin<3)
        ct = 1:size(data.Y,2);
    end

    P = zeros(net.num_neurons,num_samples);
    Z = zeros(net.num_neurons,num_samples);
    
    eta = net.etaW;
    eta_W0 = net.etaW0;
    
    for i=1:num_samples

        [net,Z_i,P(:,i)] = ...
            net.p_sample_fcn( net, data, ct(i) );
        
        Z(:,i) = Z_i;
        
        z_k = find(Z_i)';

        for k = z_k
            
            net.iteration = net.iteration+1;
            net.num_iter(k) = net.num_iter(k)+1;
            
            S_k = net.S(k,:);
            Q_k = net.Q(k,:);

            eta = net.eta*(Q_k-S_k.^2)./(exp(-S_k)+1);
            
            W_exp = exp( net.W(k,:) );
            
            net.W(k,:) = net.W(k,:) + eta.*(net.Y(:)'-W_exp)./max(eta,W_exp);

            net.S(k,:) = S_k + eta.*(net.W(k,:)-S_k);
            net.Q(k,:) = Q_k + eta.*(net.W(k,:).^2-Q_k);
        end
        
        if ~isempty(z_k)
            
            eta_w = net.eta*(net.Q0-net.S0.^2)./(exp(-net.S0)+1);

            eta_W0 = eta_w;

            net.W0 = net.W0 + eta_W0.*(Z_i-exp(net.W0))./max(eta_W0,exp(net.W0));

            net.S0 = net.S0 + eta_W0.*(net.W0-net.S0);
            net.Q0 = net.Q0 + eta_W0.*(net.W0.^2-net.Q0);
        end
    end
    
    net.etaW = eta;
    net.etaW0 = eta_W0;
end
