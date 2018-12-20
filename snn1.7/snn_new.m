function net = snn_new( num_neurons, num_inputs, varargin )
% snn_new: creates a new SNN network
%
% net = snn_new( num_neurons, num_inputs, ... )
%
% Creates a new SNN network. Optional parameters will
% modifie the network behaviour. Parameters that can
% be passed here are the same as to the
% <a href = "matlab:help snn_set">snn_set</a> function.
%
% input
%   num_neurons:  number of wta neurons.
%   num_inputs:   number of network inputs.
%
% output
%   net:          The snn network structure.
% 
% see also
%   <a href = "matlab:help snn_set">snn_set</a>
%   <a href = "matlab:help snn_train">snn_train</a>
%
% David Kappel
% 04.01.2010
%
    
    if (nargin<2)
        error( 'Not enought input arguments!' );
    end

    net.num_neurons = num_neurons;
    net.num_inputs = num_inputs;

    net.train_method = '*';
    net.sample_method = '*';
    net.performance_method = '*';
    
    net.iteration = 0;    
    net.train_time = 0;

    net = snn_set( net, varargin{:} );
end
