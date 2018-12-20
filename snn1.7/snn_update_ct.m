function [net, train_set] = snn_update_ct( net, train_set, varargin )
% snn_update_ct: update snn network with continuous time input
%
% [net, train_set] = snn_update_ct( net, train_set, ... )
%
% Update the network weights over all provided input sequences
% at once.
%
% input
%   net:              A SNN network or an integer number.
%                     If a number is passed a net is created with
%                     the given number of neurons.
%                     See <a href="matlab:help snn_new">snn_new</a>.
%   test_data:        Test data or file pattern to load data.
%
% optional arguments
%   collect:          A string containing a list of fields that should
%                     be collected from the network structure.
%                     The fields are collected after each block that was
%                     processed and stored in the data_out structure.
%                     The string must have the format:
%                     '[<field1>,<field2>,...,<fieldN>]'
%
% see also
%   <a href="matlab:help snn_new">snn_new</a>
%   <a href="matlab:help snn_load_data">snn_load_data</a>
%
% David Kappel 29.07.2011
%

    if (nargin<2)
        error('Not enought input arguments!');
    end

    if ~isstruct( net )
        error('Unexpected argument type for ''net''!');
    end

    [ collect, ...
      set_generator ] = snn_process_options( varargin, ...
                                             'collect', '', ...
                                             'set_generator', [] );
                         
    verbose = snn_options( 'verbose' );
    
    if isempty( verbose )
        verbose = true;
    end
    
    field_collectors = {};
    
    if ~isempty(collect)
        field_collectors = snn_parse_args(collect);
    end
    
    if isfield( net, 'p_required_fields' )
        req_fields = snn_parse_args(net.p_required_fields);
        field_collectors = { field_collectors{:}, req_fields{:} };
    end

    field_collectors = unique( field_collectors );
    
    if isfield( net, 'num_samples' )
        num_samples = net.num_samples;
    else
        num_samples = 1;
    end
    
    if isempty( set_generator )
        set_generator.fcn_generate = @( data_gen, i )( def_data_generator( data_gen, i ) );
        set_generator.data = train_set;
    end
    
    if (verbose)
        fprintf('updating network... trial:      0');
    end
    
    net.got_sample = false;
    
    num_trials = 0;
    
    while ~net.got_sample
        
        num_trials = num_trials + 1;
        
        if (verbose)
            fprintf('%c%c%c%c%c%5d',8,8,8,8,8,num_trials);
        end        
        
        [net, train_set] = snn_process_data( net, net.p_sample_fcn, ...
                                             set_generator, inf, num_samples, ...
                                             true, 0, field_collectors );

        net = net.p_train_fcn( net, train_set, 1:length(train_set) );
    end

    net.num_trials = num_trials;
    
    if (verbose)
        fprintf( '\n' );
    end
    
    %% dafault data generator
    function data_set = def_data_generator( data_gen, i )
        data_set = data_gen.data(i);
    end
end
