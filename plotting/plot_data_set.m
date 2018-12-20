function data_set = plot_data_set( data_set, i, fig )
% Plot the data structure.
% If a string is given the data structure
% is loaded from the file it points to.
%

    if (nargin<3)
        fig = 1;
    end

    if ischar( data_set )
        
        src_file_name = locate_data_file( data_set, i );

        fprintf( 'loading dataset %s...\n', src_file_name );
        
        if isempty( src_file_name )
            error( 'File not found!' );
        end
        
        data_set = load( src_file_name );
    end

    idx = sort_neurons_seq( data_set.sim_test, data_set.net.num_neurons );
    
    plot_data.Zt = [];
    plot_data.Xt = [];
    plot_data.Lt = [];
    plot_data.At = [];
    plot_data.Rt = [];
    plot_data.A_v = [];
    plot_data.A_w = [];
    
    plot_data.time = 0;
    plot_data.labels = [];
    
    for j = 1:min(5,length(data_set.sim_train))
        plot_data = append_data( plot_data, data_set.sim_train(j), [1,.7,.7] );
    end
    
    for j = 1:length(data_set.sim_test)
        
        train_set = data_set.sim_test{j};
        
        for k = 1:length(train_set);        
            plot_data = append_data( plot_data, train_set(k), [.7,1,.7] );
        end
        
        plot_data = append_data( plot_data, data_set.sim_free{j}, [1,1,.6] );
    end
    
    plot_data.idx = idx;
    plot_data.Zt_y_range = (1:data_set.net.num_neurons);
    
    if ~isempty( plot_data.A_v ) && ~isempty( plot_data.A_w )
        At = zeros(5,size(plot_data.At,2));
        At(1,:) = plot_data.At(1,:);
        At(2,:) = plot_data.A_v;
        At(3,:) = plot_data.A_w;
        At(4,:) = plot_data.At(1,:) - plot_data.A_v - plot_data.A_w;
        At(end,:) = plot_data.At(end,:);  
        plot_data.At = At;
    end
    
    if isfield( plot_data, 'R_all' )
        snn_plot( plot_data, '[ Lt[st], Zt[stx], At[pt], R_all[pt], _ ]', 'figure', fig );
    else        
        snn_plot( plot_data, '[ Lt[st], Zt[stx], At[pt], Rt[pt], _ ]', 'figure', fig );
    end
       
    subplot( 6, 1, 6 );
    
    plot( data_set.performance );
    title( 'performance' );
    xlabel( 'iteration' );
    
    drawnow;
end
