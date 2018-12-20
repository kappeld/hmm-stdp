function [idx, labels] = sort_neurons_seq( data_sets, num_neurons )
%
% Generate the neural sorting index over multiple sequences.
%

    spike_times = zeros(num_neurons,length(data_sets));
    num_spikes = zeros(num_neurons,length(data_sets));
    
    t_offset = 0;
    all_offsets = zeros( 1, length(data_sets) );
    
    for i = 1:length(data_sets)
        
        data = data_sets{i};
        
        for j = 3:length(data)
            
            % remove zero spikes
            tmp_Zt = data(j).Zt(:,data(j).Zt(1,:)>0);
            
            data_lenght = double( size(tmp_Zt,2) );
            
            data_spikes = sparse( double( tmp_Zt(1,:) ), ...
                                  double( 1:data_lenght ), ...
                                  double( tmp_Zt(2,:) + t_offset ), ...
                                  num_neurons, data_lenght );

            spike_times(:,i) = spike_times(:,i) + sum( data_spikes, 2 );
            num_spikes(:,i) = num_spikes(:,i) + sum( data_spikes > 0, 2 );
        end
        
        all_offsets(i) = t_offset;
        
        if isempty(data)
            continue;
        end
        
        t_offset = t_offset + data(1).time(end);
    end
    
    av_spike_times = zeros(num_neurons,1);    
    labels = zeros(num_neurons,1);
    
    for i = 1:num_neurons
        
        %n_i = find( num_spikes(i,:) > 0, 1, 'first' );
        [v,n_i] = max( num_spikes(i,:) );
        
        if isempty(n_i)
            av_spike_times(i) = inf;
        else
            av_spike_times(i) = spike_times(i,n_i)/num_spikes(i,n_i);
            
            data = data_sets{n_i};
            
            tm = data(1).time;
            lbls = data(1).labels;
            
            for j = 1:length( lbls )
                
                t_av_rel = av_spike_times(i) - all_offsets(n_i);
                
                if ( t_av_rel >= tm( lbls(j).start_sample ) ) && ...
                   ( t_av_rel <= tm( lbls(j).stop_sample ) )
                    labels( i ) = lbls(j).descriptor - 'A' + 1;
                end
            end
        end
    end

    [v,idx] = sort(av_spike_times);
    
    labels = labels(idx);
end
