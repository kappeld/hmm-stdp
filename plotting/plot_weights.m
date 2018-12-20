function [all_W, all_V, all_V0, eta_W, eta_V, eta_0] = plot_weights( path, num_iter )

    all_files = dir( [ path, '*.mat'] );
    
    if isempty( all_files )
        all_W = []; all_V = []; all_V0 = [];
        return;
    end

    data = load( [ path, all_files(1).name ] );

    num_neurons = data.net.num_neurons;
    num_inputs = data.net.num_inputs;

    if ( nargin < 2 )
        num_iter = length( all_files );
    end

    all_W = zeros( num_neurons*num_inputs, num_iter );
    all_V = zeros( num_neurons*num_neurons, num_iter );
    all_V0 = zeros( num_neurons, num_iter );

    eta_W = zeros( num_neurons*num_inputs, num_iter );
    eta_V = zeros( num_neurons*num_neurons, num_iter );
    eta_0 = zeros( num_neurons, num_iter );
    
    fprintf( 'loading data sets...  0%%' );

    for i=1:num_iter;

        data = load( [ path, all_files(i).name ] );
        
        fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*i/num_iter))

        all_W(:,i) = data.net.W(:);
        all_V(:,i) = data.net.V(:);
        all_V0(:,i) = data.net.V0(:);

        eta_W(:,i) = data.net.eta_W(:);
        eta_V(:,i) = data.net.eta_V(:);
        eta_0(:,i) = data.net.eta_0(:);

    end
    
    figure; plot( all_W' ); title( 'feedforward weights (W)' );
    figure; plot( all_V' ); title( 'lateral weights (V)' );
    figure; plot( all_V0' ); title( 'bias weights (V0)' );

    fprintf('%c%c%c%cdone.\n',8,8,8,8);
    
end
