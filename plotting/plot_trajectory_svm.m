function plot_trajectory_svm( m_data )
% eval performance of linear classifier on trajectories of
% multiple WTA circuits.
%


    samples = 0.100:0.005:0.150;
    sigma = 0.01;
    
    num_samples = length(samples);
    num_neurons = m_data.net.num_neurons;
    
    num_sets = 100;
    num_pat = 2;
    
    m_data.net.groups(2:(end))

    wta_groups = [ 1,            1, cumsum(m_data.net.groups(1:(end-1)))+1; ...
                   num_neurons,  m_data.net.groups  ];
        
    classes = [ zeros(num_samples*num_sets,1); ones(num_samples*num_sets,1) ];
    
    perf_all = zeros(1,size(wta_groups,2));
    
    idx = randperm(num_samples*num_sets*num_pat);
    
    for group = 1:size(wta_groups,2)
        
        neuron_ids = wta_groups(1,group):(wta_groups(1,group)+wta_groups(2,group)-1);

        P = zeros( length(neuron_ids), num_samples, num_sets, num_pat );

        fprintf( 'computing peth...      0%%' );

        [v,n_idx] = sort( m_data.I(neuron_ids) );

        for i = 1:num_sets
            for j = 1:num_pat
                P(:,:,i,j) = get_peth( m_data.test_data{j,i}.Zt, neuron_ids(n_idx), samples, sigma );
            end        
            fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*i/num_sets))
        end
        
        fprintf('%c%c%c%cdone.\n',8,8,8,8);
            
        X = reshape(P,[length(neuron_ids),num_samples*num_pat*num_sets]);

        x_svm = svmtrain( X(:,idx(1:min(length(idx),500))), classes(idx(1:min(length(idx),500))), ...
                          'Method', 'QP', 'boxconstraint', 20, 'QuadProg_Opts',  optimset('MaxIter', 1000 ) );

        t_classes = svmclassify(x_svm,X');
        
        perf_all(group) = sum(all(reshape( t_classes == classes, [num_samples,num_pat*num_sets] ),1))/(num_sets*num_pat);

        fprintf( '%3i   size: %3i   performance: %f\n', group, length(neuron_ids), perf_all(group) );
    end

%     colors = [ [1,0,0]; [0,1,0]; [0,0,1]; [.7,.7,0] ];
%     figure;
%     
%     for k=1:12
%         for l=(k+1):12
%             for m=(l+1):12
%                 
%                 fprintf( '%i %i %i \n', k, l, m );
% 
%                 X_svm_proj = [ x_svm.SupportVectors(k,:) - x_svm.SupportVectors(j,:); x_svm.SupportVectors(k,:) - x_svm.SupportVectors(m,:) ]';
% 
%                 hold off;
% 
%                 for j=1:num_pat
%                     for i=1:num_sets
%                         sample = zeros( num_neurons, num_samples );
%                         for c = 1:num_neurons
%                             sample(:,c) = x_svm.ScaleData.scaleFactor(c) * ...
%                                           (P(:,c,i,j)' +  x_svm.ScaleData.shift(c));
%                         end
%                         p_proj = sample'*X_svm_proj;
%                         plot( p_proj(1,:), p_proj(2,:), '-', 'Color', colors(j,:) );
%                         hold on;
%                         plot( p_proj(1,:), p_proj(2,:), '.', 'Color', colors(j,:) );
%                     end
%                 end
% 
%                 waitforbuttonpress;
%     
%             end
%         end
%     end
end
