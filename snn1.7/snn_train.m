function [net, performance, train_set] = snn_train( net, train_data, test_set, varargin )
% snn_train: train the network with given train and test data
%
% [net, performance, train_set] = snn_train( net, train_data, test_data, ... )
% [net, performance, train_set] = snn_train( net, data )
%
% Trains a SNN network with given data.
%
% input
%   net:              A SNN network or an integer number.
%                     If a number is passed a net is created with
%                     the given number of neurons.
%                     See <a href="matlab:help snn_new">snn_new</a>.
%   train_data:       Training data or file pattern to load data.
%   test_data:        Test data or file pattern to load data.
%   data:             A string pointing to a directory containing
%                     train and test data.
%
% optional arguments
%   num_runs:         number of training runs through data set
%                     default is 1.
%   step_size:        Number of samples that are trained before
%                     the performance is recalculated.
%   figure:           figure to plot performance to. Doesn't plot
%                     anything if 0 is passed.
%   keep_data_order:  If false the training data is trained in random
%                     permutation (default), else the order is kept.
%   auto_restore:     If an error occures during training the previous
%                     parameters are restored (default), else training
%                     continues without changes and the performance at
%                     current iteration is set to 'NAN'.
%                     true by default.
%   plot_rule:        Plot rule string for plotting data produced while
%                     training. See  snn_plot().
%   lead_in:          Number of samples that will be copied from previous
%                     train set.
%   collect:          A string containing a list of fields that should
%                     be collected from the network structure.
%                     The fields are collected after each block that was
%                     processed and stored in the data_out structure.
%                     The string must have the format:
%                     '[<field1>,<field2>,...,<fieldN>]'
%   [net parameters]: If 'net' is given by an interger number additional
%                     parameters may be given that are passed on to
%                     <a href="matlab:help snn_new">snn_new</a>.
%
% output
%   net:              The trained network.
%   performance:      A matrix containing the performance information
%                     for each timestep.
%   train_set:        The network training data sets.
%
% see also
%   <a href="matlab:help snn_new">snn_new</a>
%   <a href="matlab:help snn_load_data">snn_load_data</a>
%
% David Kappel 14.08.2010
%

%% init

    if (nargin<2)
        error('Not enought input arguments!');
    end

    if ~isstruct( net ) &&  ~isnumeric( net )
        error('Unexpected argument type for ''net''!');
    end

    if (nargin<3)
        if ischar(train_data)
            test_set = fullfile( train_data, 'test*.mat' );
            train_data = fullfile( train_data, 'train*.mat' );
        else
            if ~isstruct( train_data )
                error('Unexpected argument type for ''train_data''!');
            end
           
            test_set = [];
        end
    end
    
    [ train_step_size, ...
      test_step_size, ...
      fig, ...
      num_runs, ...
      keep_data_order, ...
      auto_restore, ...
      plot_rule, ...
      lead_in, ...
      collect, ...
      varargin ] = ...
        snn_process_options( varargin, ...
                             'train_step_size', 1000, ...
                             'test_step_size', inf, ...
                             'figure', [], ...
                             'num_runs', 1, ...
                             'keep_data_order', false, ...
                             'auto_restore', true, ...
                             'plot_rule', '', ...
                             'lead_in', 0, ...
                             'collect', '' );
                         
    verbose = snn_options( 'verbose' );
    
    if isempty( verbose )
        verbose = true;
    end
        
    if ischar(test_set)
        if (verbose)
            disp('loading test-set:');
        end
        test_set = snn_load_data( test_set );
    end

    performance = [];
    
    stop_button = [];
    text_box = [];
    
    test_collectors = {};
    train_collectors = {};
    
    if ~isempty(collect)
        if ~isempty( test_set )
            test_collectors = snn_parse_args(collect);
        else
            train_collectors = snn_parse_args(collect);
        end
    end

%% train network
    if (verbose)
        disp('training network:');
    end
    
    if isempty( fig )
        fig = figure;
    end

    if (nargout>2)
        [net, train_set] = snn_process_data( net, @(net,data,ct)train_network(net,data,ct), ...
                                             train_data, train_step_size, num_runs, ...
                                             keep_data_order, lead_in, train_collectors );
    else
        net = snn_process_data( net, @(net,data,ct)train_network(net,data,ct), ...
                                train_data, train_step_size, num_runs, ...
                                keep_data_order, lead_in, train_collectors );
    end
    
    delete( stop_button );
    delete( text_box );
    
    drawnow;
    
    if (verbose)
        fprintf('\n')
    end

    
%% local functions    
    
    function make_gui( out_text )
    % Set up the GUI environment.

        if isempty( stop_button )
            stop_button = uicontrol( 'Style', 'togglebutton', ...
                                     'String', 'stop',...
                                     'Position', [15 15 55 25] );

            text_box = uicontrol( 'Style', 'Text', ...
                                  'HorizontalAlignment', 'Left', ...
                                  'FontSize', 10, ...
                                  'BackgroundColor', get(fig,'Color'), ...
                                  'Position', [80 15 500 15]);
        end

        set( text_box, 'String', out_text );
    end


    function [net,Z,P] = train_network( net, data, ct )
    % Train the network

        if isnumeric(net)
            net = snn_new( net, size(data.X,1), varargin{:} );            
        end
    
        out_text = sprintf( '  %02.0f:%02.0f.%03.0f  iteration: %u', ...
                            floor(net.train_time/60), floor(mod(net.train_time+eps,60)), ...
                            floor(mod(1000*net.train_time,1000)), net.iteration );

        net_tmp = net;

        [ net, Z, P ] = net.p_train_fcn( net, data, ct );

        if isfield( data, 'time' )
            net.train_time = net.train_time + ...
                             data.time(ct(end)) - ...
                             data.time(ct(1));
        else
            net.train_time = net.train_time + ...
                             data.time_range(end) - ...
                             data.time_range(1);
        end

        if ~isempty( test_set )

            sim_net = snn_alloc( net, net.p_sample_allocators, true );
            
            [sim_net, test_set_result] = snn_process_data( sim_net, net.p_sample_fcn, ...
                                                           test_set, test_step_size, 1, ...
                                                           true, lead_in, test_collectors );
                                     
            perf = snn_performance( sim_net, test_set_result );

            if isnan(perf) && auto_restore
                it = net.iteration;
                net = net_tmp;
                net.iteration = it;
                if (verbose)
                    disp('  training error - restoring previous parameters!');
                end
                if ( numel(performance)>0 )
                    perf = performance(1,end);
                end
            end

            performance = [ performance, [perf;net.iteration] ];

            out_text = strcat( out_text, sprintf( '  performance is: %f',  perf(1) ) );

            if ( fig > 0 )

                fig = snn_plot( test_set_result(1), [ plot_rule, '_' ], ...
                                'figure', fig, varargin{:} );

                make_gui( out_text );

                plot(performance(2,:), performance(1,:), '-b');
                title( 'training network' );
                xlabel('iteration');
                ylabel('performance');

                drawnow;
            end
        else
            if ~isempty( plot_rule ) && ( fig > 0 )

                fig = snn_plot( data, plot_rule, 'figure', fig, varargin{:} );                        
                make_gui( out_text );
                
                drawnow;
            end
        end

        if (verbose)
            disp( out_text );
        end

        if ~isempty( stop_button ) && get( stop_button, 'Value' )
            
            net.p_user_interupt = true;
            if (verbose)
                disp( 'User interrupt...' );
            end
        end
    end
end
