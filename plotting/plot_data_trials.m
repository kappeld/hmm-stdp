function plot_data_trials( sim_data, num_neurons, num_inputs, idx, run_id, lbls, pat_labels, colors )
% Plot the data structure.
% If a string is given the data structure
% is loaded from the file it points to.
%

    if (nargin < 5)
        run_id = 1;
    end

    if (nargin < 6)
        lbls = get_neuron_labels( sim_data(run_id), num_neurons, 5 );
    end

    plot_heights = [ .30, .30 ];

    sim_data(run_id).Zt_y_range = 1:num_neurons;
    sim_data(run_id).idx = idx;
    
    
    h_fig = figure;
    hist_w = .08;
    x_l_border = .18;
    x_r_border = .10;
    y_border = .12;
    hist_border = .04;
    plot_border = .14;
    hist_max = 20;
    max_x_neurons = min(num_inputs,20);
    label_dist = 4;
    title_dist = 8;
    
    p_x = x_l_border + hist_w + hist_border;
    p_y = y_border + plot_heights(2) + plot_border;
    plot_w = 1 - hist_w - x_l_border - x_r_border - hist_border;    
    
    pos = get( h_fig, 'Position' );
    pos(3) = 280; pos(4) = 280;
    set( h_fig, 'Position', pos, 'Renderer', 'painters' );
    
    if ~isempty( lbls )
        sim_data(run_id).Zt = [ sim_data(run_id).Zt(1,:); ...
                                lbls( sim_data(run_id).Zt(1,:) ); ...
                                sim_data(run_id).Zt(2,:) ];
    end
    
    Lt = sim_data(run_id).Lt;
    Xt = sim_data(run_id).Xt;
    x_idx = Lt(1,:) <= max_x_neurons;
    Lt = [ Lt(1,x_idx); Lt(2,x_idx); Lt(3,x_idx) ];
    Xt = [ Xt(1,x_idx); Xt(2,x_idx) ];
    sim_data(run_id).Lt = Lt;
    sim_data(run_id).Xt = Xt;
    sim_data(run_id).Lt_y_range = 1:max_x_neurons;
    
    if ~isempty( pat_labels )
        sim_data(run_id).Zt_labels = get_data_labels( sim_data(run_id).Zt, pat_labels );
    end
    
    [h_fig,h_axes] = snn_plot( sim_data(run_id), '[ _, Lt[stlc], _, Zt[sx] ]', ...
                               'figure', h_fig, 'numbering', 'none', 'colors', colors, ...
                               'axis_grid', [ x_l_border, p_y, hist_w, plot_heights(1); ...
                                              p_x, p_y, plot_w, plot_heights(1); ...
                                              x_l_border, y_border, hist_w, plot_heights(2); ...
                                              p_x, y_border, plot_w, plot_heights(2) ] );
                                          
                                          
    axes(h_axes(2));
    h_title = title( 'network inputs' );
    
    set( h_title, 'Units', 'points' );
    pos = get( h_title, 'Position' );
    pos(2) = pos(2) + title_dist;
    set( h_title, 'Position', pos );
    
    set( h_title, 'Units', 'normalized' );
    pos = get( h_title, 'Position' );
    pos(1) = .5 - ((hist_w + hist_border)/2)/plot_w;
    set( h_title, 'Position', pos );
    
    set( gca, 'YTick', [], 'YTickLabel', {}, ...
              'YTickMode', 'manual', 'YTickLabelMode', 'manual' );

    axes(h_axes(4));
    h_title = title( 'network outputs' );
    
%     set( h_title, 'Units', 'points' );
%     pos = get( h_title, 'Position' );
%     pos(2) = pos(2) + title_dist;
%     set( h_title, 'Position', pos );
    
    set( h_title, 'Units', 'normalized' );
    pos = get( h_title, 'Position' );
    pos(1) = .5 - ((hist_w + hist_border)/2)/plot_w;
    set( h_title, 'Position', pos );
    
    xlabel( 'time [ms]' );
    set( gca, 'YTick', [], 'YTickLabel', {}, ...
              'YTickMode', 'manual', 'YTickLabelMode', 'manual' );

    axes(h_axes(1));
    sorted_hist_plot( sim_data(run_id), max_x_neurons, 1:max_x_neurons, 'Xt' );
    set( gca, 'FontName', snn_options( 'FontName' ), 'FontSize', snn_options( 'FontSize' ) );
    h_y_label = ylabel( 'input neuron' );
    %pos = get( h_y_label, 'Position' );
    %pos(1) = -label_dist;
    %set( h_y_label, 'Position', pos );
    xlabel( '' );
    set( gca, ... %'YTick', [0.5, num_inputs+0.5], ...
              ... %'YTickLabel', {'1',num2str(num_inputs)}, ...
              'YTickMode', 'manual', 'YTickLabelMode', 'manual', ...
              'XTick', [0, hist_max], ...
              'XTickLabel', {'0', sprintf('%i',hist_max)}, ...
              'XTickMode', 'manual', 'XTickLabelMode', 'manual' );
    xlim( [0, hist_max] );

    axes(h_axes(3));
    sorted_hist_plot( sim_data(run_id), num_neurons, idx, 'Zt' );
    set( gca, 'FontName', snn_options( 'FontName' ), 'FontSize', snn_options( 'FontSize' ) );
    h_y_label = ylabel( 'trial' );
    %pos = get( h_y_label, 'Position' );
    %pos(1) = -label_dist;
    %set( h_y_label, 'Position', pos );
    xlabel( '# spikes' );
    set( gca, ... %'YTick', [0.5, num_neurons+0.5], ...
              ... %'YTickLabel', {'1',num2str(num_neurons)}, ...
              'YTickMode', 'manual', 'YTickLabelMode', 'manual', ...
              'XTick', [0, hist_max], ...
              'XTickLabel', {'0', sprintf('%i',hist_max)}, ...
              'XTickMode', 'manual', 'XTickLabelMode', 'manual' );
    xlim( [0, hist_max] );

end
