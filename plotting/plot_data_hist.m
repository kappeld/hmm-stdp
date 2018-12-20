function plot_data_hist( data_set, varargin )
% Plot the data structure. If a string is given the data structure is
% loaded from the file it points to.
%

    [ num_neurons, ...
      idx, ...
      lbls, ...
      colors, ...
      pat_labels, ...
      z_labels, ...
      plot_axis, ...
      hist_axis, ...
      groups, ...
      run_id, ...
      seq_id, ...
      set_name, ...
      varargin ] = snn_process_options( varargin, ...
                                        'num_neurons', [], ...
                                        'neuron_order', [], ...
                                        'neuron_labels', [], ...
                                        'colors', [], ...
                                        'pat_labels', [], ...
                                        'z_labels', [], ...
                                        'plot_axis', [], ...
                                        'hist_axis', [], ...
                                        'groups', [], ...
                                        'run_id', 1, ...
                                        'seq_id', 1, ...
                                        'set_name', 'sim_test' );


    if isempty(num_neurons)
        if ~isempty( idx )
            num_neurons = length( idx );
        else
            num_neurons = data_set.net.num_neurons;
        end
    end
    
    if isempty(idx)
    	idx = sort_neurons_seq( data_set.(set_name), num_neurons );
    end
    
    if ( length(idx) < data_set.net.num_neurons )
        idx = [ idx, find(  not(sparse( 1, idx, true, 1, data_set.net.num_neurons )) ) ];
    end

    if isempty(lbls)
        lbls = get_neuron_labels( data_set.(set_name)(run_id), num_neurons, 5 );
    end
    
    if isempty(colors)
        if isfield(data_set,'colors')
            colors = data_set.colors;
        else
            colors = snn_options( 'colors' );
        end
    end

%% prepare data for plotting

    if (run_id > 1) && (size(data_set.(set_name),2) == 1)
        sim_data = data_set.(set_name){seq_id}(run_id);
    else
        sim_data = data_set.(set_name){seq_id,run_id};
    end

    sim_data.Zt_y_range = 1:num_neurons;
    sim_data.idx = idx;

    if ~isempty( lbls )
        sim_data.Zt = [ sim_data.Zt(1,:); ...
                        lbls( sim_data.Zt(1,:) ); ...
                        sim_data.Zt(2,:) ];
    end
    
    grp_lbls = struct;
    
    for i=2:2:(length(groups)-1)
        grp_lbls(i/2).show_border = false;
        grp_lbls(i/2).color = 1;
        grp_lbls(i/2).start_time = sim_data.time(1);
        grp_lbls(i/2).stop_time = sim_data.time(end)+0.01;
        grp_lbls(i/2).range = [ sum(groups(1:(i-1)))+1, sum(groups(1:i)) ];
    end
    
    if isempty( z_labels )
        sim_data.Zt_labels = get_data_labels( sim_data.Zt, pat_labels, grp_lbls );
    else
        sim_data.Zt_labels = z_labels;
    end
    
    if isempty( plot_axis )
        figure;
        plot_axis = gca;
    end
    
    snn_plot( sim_data, '[ Zt[stclx] ]', ...
              'axis', plot_axis, 'numbering', 'none', ...
              'colors', colors, ...
              'show_titles', false, varargin{:} );
                                          
                                          
    if ~isempty(hist_axis)
        
        axes(hist_axis);
        sorted_hist_plot( sim_data, num_neurons, idx, 'Zt' );
        set( gca, 'LineWidth', 0.8, 'FontName', snn_options( 'FontName' ), ...
                  'FontSize', snn_options( 'FontSize' ), ...
                  'YTickMode', 'manual', 'YTickLabelMode', 'manual', ...
                  'YTick', [0.5, num_neurons+0.5], ...
                  'YTickLabel', {'1', sprintf('%i',num_neurons)} );
    end

end
