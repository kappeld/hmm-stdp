function data_set = plot_trajectory( data_set, varargin )
% Plot the data structure.
% If a string is given the data structure
% is loaded from the file it points to.
%


    if ischar( data_set )
        data_set = load( data_set );
    end
    
    [ h, ...
      show_trace, ...
      smoothing, ...
      time_range, ...
      group_id ] = snn_process_options( varargin, ...
                                        'axes', [], ...
                                        'show_trace', false, ...
                                        'smoothing', 0.25, ...
                                        'time_range', [], ...
                                        'group_id', [] );
    
    colors = snn_options( 'colors' );

    g_idx = [0,cumsum( data_set.net.groups )];

    if isempty(h)
        figure;
    else
        axes(h);
    end
    
    for i=1:size(data_set.psps_evoced,1)
        
        if isempty(group_id)
            I_g = data_set.I;
            nneur = data_set.net.num_neurons;
        else
            I_g = data_set.I(and( data_set.I>g_idx(group_id), data_set.I<=g_idx(group_id+1) ));
            nneur = data_set.net.groups(group_id);
        end
        
        pc = [sin((0:(nneur-1))*(2*pi/(nneur-1)));cos((0:(nneur-1))*(2*pi/(nneur-1)))];
        
        for j=1:25
            % remove first few samples. They are very simmilar
            % since the initial conditions where the same and
            % would fudge the results.
            
            if isempty(time_range)
                t_r = 1:size(data_set.psps_evoced{i,j},2);
            else
                t_r = time_range;
            end
            
            tmp_psp = data_set.psps_evoced{i,j};
            data_set.psps_evoced{i,j} = tmp_psp(:,t_r);            
        end
        
        psp = [data_set.psps_evoced{i,1:25}];
        poj = pc*psp(I_g,:);
        p = ashape(poj(1,:),poj(2,:),smoothing,'-g');
        patch( p.x( p.apatch{1} ), p.y( p.apatch{1} ), colors(i,:) );        
        hold on;
        drawnow;
    end
    
    if show_trace
        
        if isempty(group_id)
            I_g = data_set.I;
            nneur = data_set.net.num_neurons;
        else
            I_g = data_set.I(and( data_set.I>g_idx(group_id), data_set.I<=g_idx(group_id+1) ));
            nneur = data_set.net.groups(group_id);
        end
        
        pc = [sin((0:(nneur-1))*(2*pi/(nneur-1)));cos((0:(nneur-1))*(2*pi/(nneur-1)))];
        psp = data_set.psps_evoced{1,1};
        poj = pc*psp(I_g,:);
        plot( poj(1,:),poj(2,:), '.', 'MarkerFaceColor', colors(i,:), 'MarkerSize', 1 );
        hold on;
    end
    
end
