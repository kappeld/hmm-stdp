function save_path = do_learning_task( exp_path, pat_sequences, varargin )
% DO_LEARNING_TASK runs an experiment that uses the HMM SEM network
% 
% function save_path = do_learning_task( exp_path, pat_sequences, ... )
%
% Train patterns into a HMM SEM network. Returns path that contains the
% training results.
%
% 02.11.2011
% David Kappel
%

%% init

    [ num_neurons, ...
      num_inputs, ...
      pattern_length, ...
      length_std, ...
      noise_length, ...
      input_rate, ...
      output_rate, ...
      num_train_sets, ...
      targets_per_pattern, ...
      num_patterns, ...
      num_samples, ...
      num_runs, ...
      num_train_samples, ...
      free_run_time, ...
      free_run_seqs, ...
      free_run_pat_lenghts, ...
      seq_probs, ...
      seq_ids, ...
      rec_delay, ...
      eta, ...
      save_interval, ...
      free_run_noise, ...
      iw_track_speed, ...
      w_rf, ...
      use_iw, ...
      reset_psps, ...
      use_variance_tracking, ...
      self_inhibition, ...
      tau_x_r, ...
      tau_z_r, ...
      tau_x_f, ...
      tau_z_f, ...
      tau_rf, ...
      train_method, ...
      sample_method, ...
      show_fig, ...
      pat_labels, ...
      pat_alpha, ...
      pat_beta, ...
      num_epochs, ...
      train_set_generator, ...
      performance_method, ...
      pattern_type, ...
      patterns, ...
      pretrain_ff_weights, ...
      field_collectors, ...
      net_V, ...
      input_process, ...
      regen_patterns, ...
      regen_seqs, ...
      groups, ...
      base_path, ...      
      changelog_flag, ...
      varargin ] = snn_process_options( varargin, ...
                                        'num_neurons', 10, ...
                                        'num_inputs', 200, ...
                                        'pattern_length', 0.075, ...
                                        'length_std', 0.000, ...
                                        'noise_length', 0, ...
                                        'input_rate', 100, ...
                                        'output_rate', 200, ...
                                        'num_train_sets', 1000, ...
                                        'targets_per_pattern', 10, ...
                                        'num_patterns', [], ...
                                        'num_samples', 10, ...
                                        'num_runs', 1, ...
                                        'num_train_samples', 3, ...
                                        'free_run_time', 0.500, ...
                                        'free_run_seqs', [], ...
                                        'free_run_pat_lenghts', [], ...
                                        'seq_probs', [], ...
                                        'seq_ids', [], ...
                                        'rec_delay', 0.000, ...
                                        'eta', [], ...
                                        'save_interval', 100, ...
                                        'free_run_noise', 0.00, ...
                                        'iw_track_speed', 0.001, ...
                                        'w_rf', 0, ...
                                        'use_iw', false, ...
                                        'reset_psps', false, ...
                                        'use_variance_tracking', false, ...
                                        'self_inhibition', 0, ...
                                        'tau_x_r', 0.002, ...
                                        'tau_z_r', 0.002, ...
                                        'tau_x_f', 0.02, ...
                                        'tau_z_f', 0.02, ...
                                        'tau_rf', 0.005, ...
                                        'train_method', 'rs', ...
                                        'sample_method', 'ctrf', ...
                                        'show_fig', false, ...
                                        'pat_labels', [], ...
                                        'pat_alpha', 0.2, ...
                                        'pat_beta', 0.8, ...
                                        'num_epochs', 1, ...
                                        'train_set_generator', [], ...                                        
                                        'performance_method', 'ct', ...
                                        'pattern_type', 'beta', ...
                                        'patterns', [], ...
                                        'pretrain_ff_weights', false, ...
                                        'collect', '[At,R]', ...
                                        'net_V', [], ...
                                        'input_process', 'poisson', ...
                                        'regen_patterns', false, ...
                                        'regen_seqs', false, ...
                                        'groups', [], ...
                                        'base_path', '', ...
                                        'changelog_flag', [] );


    save_path = gen_results_path( [base_path,exp_path], changelog_flag );
    
    if isempty( save_path )
        return;
    end

    data_set_files = dir( [ save_path, 'data_set_*.mat' ] );

    continue_training = 0;
    cur_iteration = 1;
    cur_epoch = 1;
    
    set_eta = true;
    
    if isempty(eta)
        eta = 0.001;
        set_eta = false;
    end

    if ~isempty( data_set_files )

        continue_training = sscanf( data_set_files(end).name, 'data_set_%d.mat' );

        u_input = input( sprintf( 'continue training at iteration? (default %i): ', continue_training ), 's' );

        if ~isempty( str2num( u_input ) )
            continue_training = str2num( u_input );
            cur_iteration = mod( continue_training, 100000 )+1;
            cur_epoch = floor( continue_training/100000 )+1;
            
            if (cur_iteration <= 0)
                cur_iteration = 1;
            end
            if (cur_epoch <= 0)
                cur_epoch = 1;
            end
        end
    end

%% create data sets

    if isempty( num_patterns )
        if ~isempty( pat_sequences )
            num_patterns = max( [ pat_sequences{:} ] );
        elseif ~isempty( train_set_generator )
            num_patterns = length( train_set_generator.pat_labels );
        else
            error( 'invalid parameters: unknown number of patterns!' );
        end
    end

    if ~isempty( pat_sequences ) 

        if ( continue_training == 0 )

            if isempty( patterns )
                patterns = gen_patterns();
            end
            
            performance = nan( num_epochs, ceil(num_train_sets/save_interval) );
            num_trials = nan( num_epochs, num_train_sets );
            all_seq_ids = zeros( num_epochs, num_train_sets );

        else

            old_path = locate_data_file( save_path, continue_training );

            if isempty( old_path )
                error( 'File not found!' );
            end

            fprintf( 'restoring experiment state %s...\n', old_path );
            old_data = load( old_path, 'patterns', 'pat_sequences', 'net', ...
                             'train_set_generator', 'performance', 'num_trials', 'all_seq_ids' );

            net = old_data.net;
            performance = old_data.performance;
            
            if isfield( old_data, 'train_set_generator' )
                train_set_generator = old_data.train_set_generator;
            end

            if isfield( old_data, 'patterns' )
                patterns = old_data.patterns;
            end
            
            if isfield( old_data, 'pat_sequences' )
                pat_sequences = old_data.pat_sequences;
            end
        end

        if isempty( pat_labels )
            if isfield( train_set_generator, 'pat_labels' )                
                pat_labels = train_set_generator.pat_labels;
            else
                num_pats = max( [pat_sequences{:}] );
                pat_labels = mat2cell( char( 'A' + (1:num_pats) - 1 ), 1, ones(1,num_pats) );
            end
        end

        num_seqs = length(pat_sequences);

        if isempty( train_set_generator )
            train_set_generator = struct();
            train_set_generator.patterns = patterns;
            train_set_generator.pat_sequences = pat_sequences;
            train_set_generator.pattern_length = pattern_length;
            train_set_generator.targets_per_pattern = targets_per_pattern;
            train_set_generator.length_std = length_std;
            train_set_generator.fcn_generate = @( data_gen, i )( generate_pattern_sequence( data_gen, i ) );
            train_set_generator.process = input_process;
        end

        if ischar( pat_sequences ) && strcmp( pat_sequences, 'basic_set_generator' )
            if isempty( train_set_generator )
                error( 'Basic train set generator required!' )
            end

            train_set_generator.patterns = patterns;
            train_set_generator.pattern_length = pattern_length;
            train_set_generator.targets_per_pattern = targets_per_pattern;
            train_set_generator.length_std = length_std;
            train_set_generator.process = input_process;
        end
        
    elseif isempty( train_set_generator )
        error( 'Argument ''train_set_genrator'' or ''pat_sequences'' required!' );
    elseif isfield( train_set_generator, 'pat_sequences' )
        num_seqs = length(train_set_generator.pat_sequences);
    else
        num_seqs = num_train_samples;
    end
    
    if isempty( seq_probs )
        seq_probs = ones( num_seqs, 1 )./num_seqs;
    end

    if ~isfield( train_set_generator, 'pat_labels' )
        train_set_generator.pat_labels = pat_labels;
    end
    
    if ~isempty( free_run_pat_lenghts )
        free_run_pat_lenghts(end+1) = free_run_time - sum(free_run_pat_lenghts);
        free_set_generator.pattern_length = free_run_pat_lenghts;
    end


%% create network

    snn_options( 'verbose', true );

    for epoch = cur_epoch:num_epochs

        if ( continue_training == 0 )

            net = snn_new( num_neurons, num_inputs, ...
                           'train_method', train_method, ...
                           'sample_method', sample_method, ...
                           'performance_method', performance_method, ...
                           'lambda', output_rate, ...
                           'self_inhibition', self_inhibition, ...
                           'iw_track_speed', iw_track_speed, ...
                           'tau_x_r', tau_x_r, ...
                           'tau_z_r', tau_z_r, ...
                           'tau_x_f', tau_x_f, ...
                           'tau_z_f', tau_z_f, ...
                           'tau_rf', tau_rf, ...
                           'w_rf', w_rf, ...
                           'mean_rec_delay', rec_delay, ...                       
                           'use_iw', use_iw, ...
                           'eta', eta, ...
                           'groups', groups, ...
                           'use_variance_tracking', use_variance_tracking, ...
                           'num_samples', num_samples, ...
                           varargin{:} );

            if regen_patterns
                train_set_generator.patterns = gen_patterns();
            end
            
            % create training/testing sets
            if regen_seqs || (epoch == 1)
                fprintf( 'creating test sets...  0%%' );

                num_test_sets = min( num_seqs, 5 );

                test_set = cell(1,num_test_sets);

                for i = 1:num_test_sets

                    test_set{i}.seq_id = i;        
                    train_set_generator.seq_id = i;

                    test_set{i}.data = train_set_generator.fcn_generate( train_set_generator, 1 );

                    fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*i/num_seqs))
                end

                fprintf('%c%c%c%cdone.\n',8,8,8,8);
                fprintf( 'creating free run sets...  0%%' );

                free_run_data = cell(1,num_test_sets);

                free_set_generator = train_set_generator;

                if ~isempty( free_run_seqs )
                    for i=1:length(free_run_seqs)
                        free_run_seqs{i}(end) = length( free_set_generator.patterns );
                    end
                    free_set_generator.free_run_seqs = free_run_seqs;
                    free_set_generator.pat_sequences = free_run_seqs;
                    free_set_generator.process = input_process;
                end

                for i = 1:num_test_sets
                    free_run_data{i}.seq_id = i;
                    free_set_generator.seq_id = i;
                    free_run_data{i}.data = free_set_generator.fcn_generate( free_set_generator, 1 );
                    fprintf('%c%c%c%c%3d%%',8,8,8,8,round(100*i/num_seqs))
                end

                fprintf('%c%c%c%cdone.\n',8,8,8,8);
            end
                       
            if pretrain_ff_weights
                pats = [ train_set_generator.patterns{1:num_patterns} ];
                net.W = repmat( log(pats/10), 1, ceil(num_neurons/num_patterns) )';
                net.W = max( -50, net.W(1:num_neurons,1:num_inputs) );
            end
            
            if ~isempty( net_V )
                net.V = net_V;
            end
        end

        continue_training = 0;

        net.num_runs = 1;
        net.hX_init = zeros( num_inputs, 2 );
        net.hZ_init = zeros( num_neurons, 2 );
        net.reset_psps = reset_psps;
        
        if set_eta
            net.eta = eta;
        end

%% train network

        tic;
        
        sim_train = {};

        for i = cur_iteration:num_train_sets

            iteration = net.iteration;
            
            if ( mod(iteration,save_interval) == 0 )
                
                net_sim = net;

                for j = 1:num_test_sets

                    net_sim.use_inhibition = true;
                    [sim_test{j},net_sim] = snn_simulate( net_sim, test_set{j}.data, 'collect', field_collectors );
                    [sim_free{j},net_sim] = snn_simulate( net_sim, free_run_data{j}.data, 'collect', field_collectors );        
                end

                performance( epoch, iteration/save_interval+1 ) = snn_performance( net, sim_test{1} );
                fprintf( 'performance is: %f\n', performance(iteration/save_interval+1) );

                toc; tic;

                data_set.net = net;
                data_set.sim_train = sim_train;
                data_set.sim_test = sim_test;
                data_set.sim_free = sim_free;
                data_set.performance = performance;

                if show_fig
                    plot_data_set( data_set );
                end

                file_name = [ save_path, sprintf( 'data_set_%03d%05d.mat', epoch-1, iteration ) ];

                fprintf( 'saving results to: %s...\n', file_name );

                save( file_name, 'net', 'sim_train', 'sim_test', 'sim_free', 'all_seq_ids', ...
                      'performance', 'train_set_generator', 'free_set_generator', ...
                      'num_trials', 'seq_probs' );
            end
            
            iteration = iteration+1;

            if isempty( seq_ids )
                seq_id = find( cumsum(seq_probs) > rand(), 1 );
            else
                seq_id = seq_ids( mod(i-1,length(seq_ids))+1 );
            end

            all_seq_ids( epoch, i ) = seq_id;
            train_set_generator.seq_id = seq_id;

            if isfield( train_set_generator, 'pat_sequences' )
                fprintf( '\niteration: %i, sequence %s\n', iteration, ...
                         [ pat_labels{ train_set_generator.pat_sequences{seq_id} } ] );
            else
                fprintf( '\niteration: %i\n', iteration );
            end
            
            net.num_runs = num_runs;

            [net,sim_train] = snn_update_ct( net, [], ...
                                             'set_generator', train_set_generator, ...
                                             'collect', field_collectors );

            num_trials( epoch, i ) = net.num_trials;

            net.num_runs = 1;

            sim_test = cell(num_seqs,1);
            sim_free = cell(num_seqs,1);

            fprintf( 'R_all: [ ' );
            fprintf( '%f ', sim_train(:).R );
            fprintf( '] / [ ' );
            fprintf( '%f ', wta_softmax( [sim_train(:).R]' ) );
            fprintf( '] R_mean %f\n', net.r_mean );
        end

        cur_iteration = 1;
        
    end
    
    function pats_ = gen_patterns()
    % generate input patterns
        
        pats_ = cell(1,num_patterns+1);

        for i_ = 1:num_patterns
            switch pattern_type
                case 'beta'
                    pat_beta_dist_mean = (pat_alpha/(pat_alpha+pat_beta));
                    max_in_rate = input_rate/pat_beta_dist_mean;
                    pats_{i_} = max_in_rate*betarnd( pat_alpha, pat_beta, [num_inputs, 1] );
                case 'beta_normalised'
                    pats_{i_} = input_rate*num_inputs*...
                        mk_normalised( betarnd( pat_alpha, pat_beta, [num_inputs, 1] ) );
                case 'sparse'
                    pat = zeros( num_inputs, 1 );
                    pat( (i_):num_patterns:num_inputs ) = input_rate;
                    pats_{i_} = pat;
                otherwise
                    error( 'Unknown pattern type!' );
            end
        end

        pats_{ num_patterns + 1 } = free_run_noise * input_rate * ones( num_inputs, 1 );
    end
end
