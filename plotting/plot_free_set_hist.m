function data_set = plot_free_set_hist( data_set, seq_id, lbls )
% Plot the data structure.
% If a string is given the data structure
% is loaded from the file it points to.
%

    idx = sort_neurons_seq( data_set.sim_test, data_set.net.num_neurons );

    if isfield( data_set, 'train_set_generator' )
        data_set.pat_labels = data_set.train_set_generator.pat_labels;
    end
    
    % FIXMEEE
    data_set.pat_labels = { 'A', 'B', 'C', 'D', 'E', 'a', 'b', 'c', 'd', '' };
    
    if (nargin < 3)
        lbls = get_neuron_labels( data_set.sim_test, data_set.net.num_neurons, data_set.pat_labels );
    end
    
    sim_data = data_set.sim_free{seq_id};
    
    for i = length( sim_data.labels ):-1:1
        
        if ( sim_data.labels(i).stop_sample - sim_data.labels(i).start_sample ) > 1
            break;
        end
    end
    
    sim_data.labels = sim_data.labels(1:i);

    if isfield( data_set, 'colors' )
        colors = data_set.colors;
    else
        colors =   [ 0.0, 0.0, 1.0; ... %A
                     0.0, 0.7, 0.0; ... %B
                     1.0, 0.0, 0.0; ... %C
                     0.9, 0.9, 0.0; ... %D
                     0.7, 0.7, 0.7; ... %E
                     0.8, 0.0, 0.8; ... %F
                     0.0, 0.6, 0.8; ... %G
                     0.8, 0.2, 0.2; ... 
                     0.2, 0.8, 0.2; ...
                     0.0, 0.0, 0.0; ];
    end
    
    if isfield( data_set.net, 'groups' )
        groups = data_set.net.groups;
    else
        groups = struct;
    end
    
    plot_data_hist( sim_data, data_set.net.num_neurons, ...
                    data_set.net.num_inputs, idx, 1, lbls, ...
                    data_set.pat_labels, colors, groups );
    
    drawnow;
end
