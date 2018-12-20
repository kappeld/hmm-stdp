function plot_trajectory_jpca( m_data )
% Plot the data structure.
% If a string is given the data structure
% is loaded from the file it points to.
%

    % params:
    %neurons = 126:148;
    neurons = 68:80;
    %neurons = 1:m_data.net.num_neurons;

    samples = 0.105:0.001:0.145;
    sigma = 0.01;
    subset_size = 10;
    
    [num_pat,num_sets] = size(m_data.test_data);
    num_samples = length(samples);
    num_neurons = length(neurons);
    
    num_sets = 20;
        
    P = zeros( num_neurons, num_samples, num_pat, num_sets );
    
    pca_coeff = [sin((0:(num_neurons-1))*(2*pi/(num_neurons-1))); ...
                 cos((0:(num_neurons-1))*(2*pi/(num_neurons-1)))]';
    
    fprintf( 'computing peth...      0%%' );
    
    [v,idx] = sort( m_data.I(neurons) );
    
    for i = 1:num_sets
        for j = 1:num_pat
            P(:,:,j,i) = get_peth( m_data.test_data{j,i}.Zt, neurons(idx), samples, sigma );
        end        
        fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*i/num_sets))
    end
    
    
    fprintf('%c%c%c%cdone.\n',8,8,8,8);
    
    colors = [ [1,0,0]; [0,1,0]; [0,0,1]; [.7,.7,0] ];

    P_jpca = jPCA( reshape(P,[num_neurons,num_samples,num_pat*num_sets]) );
    
    figure; hold on;
    
    for j=1:num_pat
        for i=0:(subset_size-1)
            plot( P_jpca(1,:,i*num_pat+j), P_jpca(2,:,i*num_pat+j), '-', 'Color', colors(j,:) );
            plot( P_jpca(1,1,i*num_pat+j), P_jpca(2,1,i*num_pat+j), '.', 'Color', colors(j,:) );
            plot( P_jpca(1,end,i*num_pat+j), P_jpca(2,end,i*num_pat+j), 'r', 'Color', colors(j,:) );
        end
    end
end
