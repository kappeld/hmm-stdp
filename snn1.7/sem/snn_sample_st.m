function [net,Z,P] = snn_sample_st( net, data, ct )
% Samples a spike from sem network in every time step.
%
% [net,Z,P] = wta_sample_st( net, data, ct )
%
% Generates the output of a given wta-network by drawing
% a winner neuron from the current output probabilities
% (standard version).
%
% inputs:
%   net:  A wta-network, see: wta-new()
%   data: A data structure to be simulated.
%   ct:   Current time index.
%
% output:
%   net:      The (modified) network structure
%   Z:        Network output spikes
%
% @fields:
%   U    zeros    [num_neurons,1]    Membrane potentials
%   A    zeros    [1,1]              Network activation
%

    Y=data.X(:,ct);
    U = net.W*Y; 
    [P,A] = wta_softmax(U);
    Z = wta_draw( P );
    net.Y = Y(:,end);
    net.U = U(:,end);
    net.A = A(:,end);
end
