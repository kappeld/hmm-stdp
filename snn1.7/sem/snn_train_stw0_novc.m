function [net,Z,P] = snn_train_stw0_novc( net, data, ct )
% Stochastic EM learning with variance tracking and bias weight update.
%
% [net,Z,P] = snn_train_st_novc( net, data, ct )
%
%
% @parameters:
%   eta                 0.01    learn rate
%   sample_method       'st'    default sample method
%   performance_method  ce      default performance method
%   temperature         -0.01   initial weight temperature
%
% @fields:
%   W         rand[temperature]  [num_neurons,num_inputs]  network weights
%   W0        rand[temperature]  [num_neurons,1]           network bias
%   num_iter  zeros              [num_neurons,num_inputs]  # of iterations
%

    num_samples = length(ct);

    if (nargin<3)
        ct = 1:size(data.Y,2);
    end

    P = zeros(net.num_neurons,num_samples);
    Z = zeros(net.num_neurons,num_samples);
    
    for i=1:num_samples

        [net,Z_i,P(:,i)] = ...
            net.p_sample_fcn( net, data, ct(i) );
        
        Z(:,i) = Z_i;
        
        z_k = find(Z_i)';

        for k = z_k
            
            net.iteration = net.iteration+1;
            net.num_iter(k) = net.num_iter(k)+1;
            
            W_exp = exp( net.W(k,:) );
            
            net.W(k,:) = net.W(k,:) + net.eta.*(net.Y(:)'-W_exp)./max(net.eta,W_exp);
        end
        
        if ~isempty(z_k)
            
            net.W0 = net.W0 + net.eta.*(Z_i-exp(net.W0))./max(net.eta,exp(net.W0));
        end
    end
end
