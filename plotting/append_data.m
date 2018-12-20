function [data_set] = append_data( data_set, data, color )
% append data to the end of the dataset.

    t_offset = data_set.time(end);
    
    if isempty( data )
        return;
    end

    Zt = data.Zt; if ~isempty(Zt), Zt(2,:) = Zt(2,:) + t_offset; end;
    Xt = data.Xt; if ~isempty(Xt), Xt(2,:) = Xt(2,:) + t_offset; end;
    At = data.At; if ~isempty(At), At(end,:) = At(end,:) + t_offset; end;
    
    time = data.time(2:end) + t_offset;

    if isfield( data, 'Lt' )
        Lt = data.Lt; if ~isempty(Lt), Lt(3,:) = Lt(3,:) + t_offset; end;
    end
    
    if isfield( data, 'R' )
        cur_R = data.R; %/size( data.Zt, 2 );
        Rt = [ cur_R, cur_R; time(1), time(end) ];
    else
        if isfield( data_set, 'Rt' ) && ~isempty( data_set.Rt )
            Rt_size = size( data_set.Rt, 1 ) - 1;
        else
            Rt_size = 1;
        end
        Rt = [ nan(Rt_size,2); [time(1), time(end)] ];
    end
    
    labels.start_time = t_offset;
    labels.stop_time = time(end);
    
    if ( nargin > 2 )
        labels.color = color;
    end

    data_set.Zt = [data_set.Zt, Zt];
    data_set.Xt = [data_set.Xt, Xt];
    data_set.At = [data_set.At, At];
    data_set.Rt = [data_set.Rt, Rt];
    data_set.time = [data_set.time(1:end-1), time];
    data_set.labels = [data_set.labels, labels];
    
    if isfield( data, 'A_v' ) && isfield( data, 'A_w' )
        data_set.A_v = [data_set.A_v, data.A_v];
        data_set.A_w = [data_set.A_w, data.A_w];
    end    
    
    if isfield( data, 'Lt' )
        data_set.Lt = [data_set.Lt, Lt];
    end
    
    if isfield( data, 'R_all' )
        if ~isfield( data_set, 'R_all' )
            data_set.R_all = [];
        end
        
        data_set.R_all = [data_set.R_all, data.R_all];
    end
    
end
