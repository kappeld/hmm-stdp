function plot_free_run( data_set, seq_id )
% Plot free run data from given simulation data.
%

    if isfield( data_set, 'pat_labels' )
        pat_labels = data_set.pat_labels;
    elseif isfield( data_set.train_set_generator, 'pat_labels' )
        pat_labels = data_set.train_set_generator.pat_labels;
    else
        pat_labels = { 'A', 'B', 'C', 'D', 'E', '' };
    end

    idx = sort_neurons_seq( data_set.sim_test, data_set.net.num_neurons );

    lbls = get_neuron_labels( data_set.sim_test, data_set.net.num_neurons, pat_labels );
    
    sim_data = data_set.sim_free{seq_id};
    sim_data.idx = idx;
    sim_data.Zt_y_range = 1:data_set.net.num_neurons;

    start_time = 0;
    
    Zt_l = [ sim_data.Zt(1,:); ...
             lbls( sim_data.Zt(1,:) ); ...
             sim_data.Zt(2,:) ];
    
    sim_data.Zt = [ sim_data.Zt(1,:); ...
                    lbls( sim_data.Zt(1,:) ); ...
                    sim_data.Zt(2,:) ];

    sim_data.Zt_labels = get_data_labels( Zt_l, pat_labels );                
                
    h_fig = figure;
    
    pos = get( h_fig, 'Position' );
    pos(3) = 400; pos(4) = 120;
    set( h_fig, 'Position', pos, 'Renderer', 'painters' );

    snn_plot( sim_data, '[ Zt[stl] ]', 'numbering', 'none', 'figure', h_fig, ...
              'start_time', start_time, 'reset_start_time', true, ...
              'axis_grid', [ 0.08, 0.28, 0.87, 0.52 ] );
          
    %xlim( [0, 1000*sim_data.time(end)-start_time-100] );
          
    xlabel( 'time [ms]' );
    ylabel( 'output neuron' );
    
    drawnow;
end
