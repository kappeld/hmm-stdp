function data_set = snn_assure_continuous_time( data_set )
% assure continuous time data set.
%
% data_set = snn_assure_continuous_time( data_set )
%
% Make shure that the dataset contains at least valid fields
% Xt and Tt.
%

    if isfield( data_set, 'Xt' ) && isfield( data_set, 'Tt' )
        return;
    end

    if ( isfield( data_set, 'sample_rate' ) )
        sample_rate = data_set.samples_rate;
    elseif ~isempty( snn_options( 'sample_rate' ) )
        sample_rate = snn_options( 'sample_rate' );
    else
        sample_rate = 1000;            
    end

    if ( ~isfield( data_set, 'Xt' ) )
        if ~isfield( data_set, 'X' )
            error( 'A valid dataset must contain field ''X'' or ''Xt''!' );
        end
        
        [i,j,s] = find(data_set.X);
        data_set.Xt = [i,j/sample_rate]';
        data_set.X_range = [1, size(data_set.X,1)];
    end
    
    if ~isfield( data_set, 'Tt' ) && isfield( data_set, 'T' )
        
        [i,j,s] = find(data_set.T);
        data_set.Tt = [i,s,j/sample_rate]';
        data_set.T_range = [1, size(data_set.T,1)];
    end    
end
