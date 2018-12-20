function [data_out,net] = snn_simulate( net, data_in, varargin )
% snn_simulate: generates the output of a given network
%
% data = snn_simulate( net, data, ... )
%
% Generates the output of a given network.
%
% input
%   net:           A SNN network.
%                  See <a href="matlab:help snn_new">snn_new</a>.
%   data:          Data set or file pattern to load data.
%
% optional arguments
%   step_size:        Number of samples that are trained before
%                     the performance is recalculated.
%   lead_in:          Number of samples that will be copied from previous
%                     train set.
%   collect:          A string containing a list of fields that should
%                     be collected from the network structure.
%                     The fields are collected after each block that was
%                     processed and stored in the data_out structure.
%                     The string must have the format:
%                     '[<field1>,<field2>,...,<fieldN>]'
%
% output
%   data:          The network simulation data.
%
% see also
%   <a href="matlab:help snn_new">snn_new</a>
%   <a href="matlab:help snn_load_data">snn_load_data</a>
%   <a href="matlab:help snn_train">snn_train</a>
%
% David Kappel 12.07.2010
%

    if (nargin<2)
        error('Not enought input arguments!');
    end

    if ~isstruct( net )
        error('Unexpected argument type for ''data''!');
    end

    [ step_size, ...
      lead_in, ...
      collect, ...
      set_generator, ...
      verbose ] = ...
        snn_process_options( varargin, ...
                             'step_size', inf, ...
                             'lead_in', 0, ...
                             'collect', '', ...
                             'set_generator', [], ...
                             'verbose', snn_options( 'verbose' ) );
    
    if isempty( verbose )
        verbose = true;
    end
    
    if isempty( set_generator )
        set_generator.fcn_generate = @( data_gen, i )( def_data_generator( data_gen, i ) );
        set_generator.data = data_in;
    end
                        
    if (verbose)
        disp('simulating network...');
    end

    [net, data_out] = snn_process_data( net, net.p_sample_fcn, ...
                                        set_generator, step_size, 1, ...
                                        true, lead_in, snn_parse_args(collect) );

    %% dafault data generator
    function data_set = def_data_generator( data_gen, i )
        data_set = data_gen.data(i);
    end
end
