function [h,fig] = get_grid_layout( widths, heights, varargin )
% Crate a set of axis fitting the grid
%
% h = get_layout( widths, heights )
%
% David Kappel
% 26.01.2011
%
%

    [ fig, ...
      left_border, ...
      right_border, ...
      top_border, ...
      bottom_border, ...
      plot_x_border, ...
      plot_y_border, ...
      dx, ...
      dy ] = ...
        snn_process_options( varargin, ...
                             'figure', [], ...
                             'left_border', 32, ...
                             'right_border', 8, ...
                             'top_border', 24, ...
                             'bottom_border', 32, ...
                             'plot_x_border', 6, ...
                             'plot_y_border', 6, ...
                             'dx', 0.6, ...
                             'dy', 2 );


    if isempty(fig)
        fig = figure;
    end
    
    set( fig, 'Units', 'points' );
    
    if numel(plot_x_border) == 1
        plot_x_border = repmat(plot_x_border,1,(length(widths)-1));
    end

    if numel(plot_y_border) == 1
        plot_y_border = repmat(plot_y_border,1,(length(heights)-1));
    end
    
    plot_x_border = [plot_x_border(:);0];
    plot_y_border = [plot_y_border(:);0];
    
    plot_width = sum( widths )*dx + left_border + right_border + sum(plot_x_border);
    plot_height = sum( heights )*dy + top_border + bottom_border + sum(plot_y_border);
    
    pos = get( fig, 'Position' );    
    pos(3) = plot_width; pos(4) = plot_height;
    set( fig, 'Position', pos, 'Renderer', 'painters', 'Resize', 'off' );

    axis_grid = zeros( length(widths)*length(heights), 4 );
    
    i = 1;
    
    py = top_border;
    
    for j=1:length(heights)
        
        px = left_border;
        height = heights(j)*dy;
        
        for k=1:length(widths)
            
            width = widths(k)*dx;
            
            axis_grid(i,:) = [ px/plot_width, (plot_height-py-height)/plot_height, ...
                               width/plot_width, height/plot_height ];
            
            i = i+1;            
            px = px + width + plot_x_border(k);
        end
        
        py = py + height + plot_y_border(j);
    end

    [fig,h] = snn_plot( [], 'init', 'figure', fig, 'axis_grid', axis_grid );

end
