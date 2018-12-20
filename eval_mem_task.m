function eval_mem_task( src_path, run_id, varargin )
% run_stoch_pattern_test( src_path )
%
% test stochastic replay of patterns
%
% David Kappel
% 24.01.2012
% 

%% create data sets

    [ num_trials, ...
      dt, ...
      sigma, ...
      test_data_generator, ...
      free_data_generator, ...
      eval_second_half_only, ...
      r_min ] = snn_process_options( varargin, ...
                                        'num_trials', 100, ...
                                        'dt', 0.001, ...
                                        'sigma', 0.01, ...
                                        'test_data_generator', [], ...
                                        'free_data_generator', [], ...
                                        'eval_second_half_only', false, ...
                                        'r_min', [] );


                                    
    save_file = [ src_path, 'mem_task_test.mat' ];

    if nargin < 2
        src_file_name = locate_data_file( src_path, 'last' );
    else
        src_file_name = locate_data_file( src_path, run_id );
    end
    
    data = load( src_file_name );
    
    if ~isempty(test_data_generator)
        data_generator = test_data_generator;
    elseif isfield( data, 'train_set_generator' )
        data_generator = data.train_set_generator;
    else
        error( 'Data field ''train_set_generator'' required!' );
    end

    num_seqs = length( data_generator.pat_sequences );

    test_data = cell(num_seqs,num_trials);
    free_data = cell(num_seqs,num_trials);

    old_verbose_option = snn_options( 'verbose', false );

    net = snn_set( data.net );
    
    st = data.sim_test{1}.time(1);
    samples = st:dt:data.sim_test{1}.time(end);
    
    fil = inline('exp( -((x-mu).^2)./(2*sigma^2) )', 'x', 'mu', 'sigma');
    
    pat_labels = { data_generator.pat_labels{:}, '' };
    
%% simulate test set
    
    peth_test = cell(num_seqs,1);
    
    test_data_generator = data_generator;
    test_data_generator.time_padding = 0.050;
    
    test_data_generator.seq_id = num_seqs;
    [tmp,net] = snn_simulate( net, [], 'set_generator', test_data_generator );
    
    num_spikes_test = zeros(net.num_neurons,1);
    
    for seq_id = 1:num_seqs
        
        fprintf( 'simulating test sequence %s...      0%%', ...
                 [ pat_labels{data_generator.pat_sequences{seq_id} } ] );
             
        test_data_generator.seq_id = seq_id;
        
        cur_peth = zeros( net.num_neurons, length(samples) );
    
        for trial = 1:num_trials
            
            fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*trial/num_trials))
            
            test_data{seq_id,trial} = snn_simulate( net, [], 'set_generator', test_data_generator );

            Zt = test_data{seq_id,trial}.Zt;

            for i = 1:size( Zt, 2 )
                cur_peth(Zt(1,i),:) = cur_peth(Zt(1,i),:) + fil( samples, Zt(2,i), sigma );
                num_spikes_test(Zt(1,i)) = num_spikes_test(Zt(1,i),:) + 1;
            end
        end
        
        peth_test{seq_id} = cur_peth;
    
        fprintf('%c%c%c%cdone.\n',8,8,8,8);       
        
    end
    
%% simulate free set

    samples_free = st:dt:data.sim_free{1}.time(end);
    
    peth_free = cell(num_seqs,1);
    
    if ~isempty(free_data_generator)
        data_generator = free_data_generator;
    elseif isfield( data, 'free_set_generator' )
        data_generator = data.free_set_generator;
    else
        for seq_id = 1:num_seqs
            seq = data_generator.pat_sequences{seq_id};
            seq = [ seq(1:ceil(length(seq)/2)), length(data_generator.patterns) ];
            data_generator.pat_sequences{seq_id} = seq;
        end
        
        pattern_length = data_generator.pattern_length;
        
        if ( numel( pattern_length ) == 1 )
            pattern_length = repmat( pattern_length, length( seq ), 1 );
        else
            pattern_length = pattern_length( 1:ceil(length(seq)/2) );
        end
        
        pattern_length(end-1) = 2*pattern_length(end-1);
        pattern_length(end) = (length(seq)-2)*pattern_length(end);
        
        data_generator.pattern_length = pattern_length;
    end
    
    num_spikes_free = zeros(net.num_neurons,1);
    
    
    
    for seq_id = 1:num_seqs
        
        fprintf( 'simulating free run sequence %s...      0%%', ...
                 [ pat_labels{data_generator.pat_sequences{seq_id} } ] );
             
        data_generator.seq_id = seq_id;
        
        cur_peth = zeros( net.num_neurons, length(samples_free) );
        
        for trial = 1:num_trials

            fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*trial/num_trials))
            
            free_data{seq_id,trial} = snn_simulate( net, [], 'set_generator', data_generator );

            Zt = free_data{seq_id,trial}.Zt;
            
            num_spikes = zeros(1,net.num_neurons);
            av_time = zeros(1,net.num_neurons);

            for i = 1:size( Zt, 2 )
                cur_peth(Zt(1,i),:) = cur_peth(Zt(1,i),:) + fil( samples_free, Zt(2,i), sigma );
                num_spikes(Zt(1,i)) = num_spikes(Zt(1,i)) + 1;
                av_time(Zt(1,i)) = av_time(Zt(1,i)) + Zt(2,i);
                num_spikes_free(Zt(1,i)) = num_spikes_free(Zt(1,i),:) + 1;
            end
            
            free_data{seq_id,trial}.num_spikes = num_spikes;
            av_time(num_spikes>0) = av_time(num_spikes>0)./num_spikes(num_spikes>0);
            free_data{seq_id,trial}.av_time = av_time;
        end
        
        peth_free{seq_id} = cur_peth;
    
        fprintf('%c%c%c%cdone.\n',8,8,8,8);
        
    end
    
    snn_options( 'verbose', old_verbose_option );

    snn_options( 'verbose', old_verbose_option );
    
%% plot data

    peth_test_all = [peth_test{:}];

    [v_test,idx] = max(peth_test_all');
    idx( v_test < 1 ) = inf;
    [Y,I] = sort(idx);

    spear_free = zeros( num_seqs, num_trials );
    spear_av = zeros(1,num_seqs);
    
    rank_test = cell(1,num_seqs);
    rank_free = cell(1,num_seqs);
    
    max_firing_times = zeros( num_seqs, net.num_neurons );
    max_firing_times_free = zeros( num_seqs, net.num_neurons );
    
    for seq_id = 1:num_seqs

        [v_test,idx] = max(peth_test{seq_id}');
        max_firing_times(seq_id,:) = idx*dt;
        [Z,J] = sort(idx);
        
        rank_test{seq_id} = zeros(1,net.num_neurons);        
        rank_test{seq_id}(J) = 1:net.num_neurons;
        
        [v_free,idx] = max(peth_free{seq_id}');
        max_firing_times_free(seq_id,:) = idx*dt;
        [Z,J] = sort(idx);
        
        rank_free{seq_id} = zeros(1,net.num_neurons);        
        rank_free{seq_id}(J) = 1:net.num_neurons;
        
        total_num_spikes = zeros(1,net.num_neurons);
        total_av_time = zeros(1,net.num_neurons);
        
        for i=1:num_trials
            total_num_spikes = total_num_spikes+free_data{seq_id,i}.num_spikes;
            total_av_time = total_av_time+free_data{seq_id,i}.av_time;
        end;

        if isempty(r_min)
            r_min = num_trials; % min number of spikes
        end
        
        if eval_second_half_only
            total_num_spikes( max_firing_times_free(seq_id,:) < peth_begin_time ) = 0;
            total_num_spikes( max_firing_times_free(seq_id,:) > (data.sim_free{1}.time(end) - 0.050) ) = 0;
        end

        spear_av(seq_id) = corr( max_firing_times_free(seq_id,J(total_num_spikes(J) > r_min))', ...
                                 max_firing_times(seq_id,J(total_num_spikes(J) > r_min))', 'Type','Spearman' );
        
        for trial = 1:num_trials
            
            num_spikes = free_data{seq_id,trial}.num_spikes;
            av_time = free_data{seq_id,trial}.av_time;

            [Z,J] = sort(av_time);
            
            spear_free( seq_id, trial ) = corr( av_time( J(num_spikes(J) > 0) )', ...
                                                max_firing_times(seq_id,J(num_spikes(J) > 0))', 'Type','Spearman' );
        end
    end
    
    rank_free = cell(1,num_seqs);
    max_firing_times_free = zeros( num_seqs, net.num_neurons );
    
    for seq_id = 1:num_seqs

        [v_free,idx] = max(peth_free{seq_id}');
        max_firing_times_free(seq_id,:) = idx*dt;
        [Z,J] = sort(idx);
        rank_free{seq_id} = zeros(1,net.num_neurons);
        rank_free{seq_id}(J) = 1:net.num_neurons;
    end
    
    save( save_file, 'net', 'test_data', 'free_data', 'peth_test', 'peth_free', ...
                     'rank_test', 'rank_free', 'spear_free', 'spear_av', 'I', 'num_spikes_test', ...
                     'num_spikes_free', 'max_firing_times', 'max_firing_times_free' );
                 
    fprintf( 'results saved to file ''%s''.\n', save_file );

end
