function data_set = snn_assure_discrete_time( data_set )
% assure discrete time data set.
%
% data_set = snn_assure_discrete_time( data_set )
%
% Make shure that the dataset contains at least valid fields
% X, T, time.
%

    if isfield( data_set, 'X' ) && isfield( data_set, 'T' )
        return;
    end

    if ( ~isfield( data_set, 'sample_rate' ) )
        sample_rate = data_set.samples_rate;
    elseif ~isempty( snn_options( 'sample_rate' ) )
        sample_rate = snn_options( 'sample_rate' );
    else
        sample_rate = 1000;            
    end

    if ( ~isfield( data_set, 'X' ) )
        if ~isfield( data_set, 'Xt' )
            error( 'A valid dataset must contain field ''X'' or ''Xt''!' );
        end
        
        data_set.X = sparse(data_set.Xt(1,:),round(data_set.Xt(2,:)*sample_rate),1, ...
                            data_set.X_range(end)-data_set.X_range(1), ...
                            round( data_set.time(end)-data_set.time(1)) );
    end
    
    if ( ~isfield( data_set, 'T' ) )
        if ~isfield( data_set, 'Tt' )
            error( 'A valid dataset must contain field ''T'' or ''Tt''!' );
        end
        
        data_set.T = sparse(data_set.X(1,:),data_set.X(3,:),round(data_set.X(2,:))*sample_rate, ...
                            data_set.T_range(end)-data_set.T_range(1), ...
                            round( data_set.time(end)-data_set.time(1)) );
    end
end
