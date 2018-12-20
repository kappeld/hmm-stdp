function [labels, num_spikes] = get_neuron_labels( sim_data, num_neurons, pat_labels, mean_ft )
%
% Get the neuron lables from the porvided simulation data.
%

    num_patterns = length( pat_labels );

    num_spikes = zeros(num_neurons,num_patterns);
    
    if ~iscell( sim_data )
        tmp = cell(1,1);
        tmp{1} = sim_data;
        sim_data = tmp;
    end
    
    if (nargin < 4)
        
        mean_ft = zeros( length(sim_data), num_neurons );
        num_spikes = zeros( length(sim_data), num_neurons );
        
        for i_c = 1:length(sim_data)

            data_set = sim_data{i_c};

            if isempty( data_set )
                continue;
            end

            for i_s = 1:length( data_set )

                Zt = data_set(i_s).Zt;
                seq_id = data_set(i_s).seq_id;

                for t = 1:size( Zt, 2 )

                    mean_ft(seq_id,Zt(1,t)) = mean_ft(seq_id,Zt(1,t)) + Zt(2,t);
                    num_spikes(seq_id,Zt(1,t)) = num_spikes(seq_id,Zt(1,t)) + 1;
                end
            end
        end
        
        mean_ft = mean_ft./num_spikes;
    end
    
    labels = -ones(1,num_neurons);
    
    for i_c = 1:length(sim_data)
        
        data_set = sim_data{i_c};
        
        if isempty( data_set )
            continue;
        end
        
        for i_s = 1:length( data_set )
            
            lbls = data_set(i_s).labels;
            tm = data_set(i_s).time;
            seq_id = data_set(i_s).seq_id;

            for i = 1:num_neurons
                for j = 1:length( lbls )

                    l_id = strmatch( lbls(j).descriptor, pat_labels, 'exact' );

                    if ( mean_ft(seq_id,i) >= tm( lbls(j).start_sample ) ) && ...
                       ( mean_ft(seq_id,i) <= tm( lbls(j).stop_sample ) ) && ...
                       ( labels( i ) < 0 ) && ~isempty(l_id)
                         labels( i ) = l_id;
                    end
                end
            end
        end
    end
end
