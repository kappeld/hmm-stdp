function  [ d_W, d_V, d_V0 ] = snn_update( net, Zt, Xt )
% compute update for synaptic weights
%
% [ d_W, d_V, d_V0 ] = snn_update( net, Z, data.Xt )
%

    d_W = zeros( net.num_neurons, net.num_inputs );
    d_V = zeros( net.num_neurons, net.num_neurons );
    d_V0 = zeros( net.num_neurons, 1 );
    
    W_exp = exp(net.W);
    V_exp = exp(net.V);
    V0_exp = exp(net.V0);
    
    j = 1;
    
    for i=1:size(Zt,2)
    
        while (j < size(Xt,2)) && (Xt(2,j) < Zt(2,j))
            j = j+1;
        end
        
        a = sum( epx( abs( Xt(2,1:j-1) - Zt(2,j) )/net.tau_x ) );
        b = sum( epx( abs( Xt(2,1:j-1) - Zt(2,j) )/net.tau_x ) );
    end
end
