function options = default_options()
   options =    { 'base_path', '', ...                % base path for data storage                  
                  'pattern_length', 0.050, ...        % each pattern has constant length (s)
                  'free_run_time', 0.400, ...         % free run time (s)
                  'input_rate', 15, ...               % mean input rate per neuron (Hz)
                  'output_rate', 200, ...             % circuit output rate (Hz)
                  'sample_method', 'ctrf', ...        % continuous time spiking
                  'train_method', 'rs', ...           % continuous time batch update
                  'pattern_type', 'beta', ...         % draw pattern from beta distribution
                  'pat_alpha', 0.2, ...               % alpha parameter of beta distribution
                  'pat_beta', 0.8, ...                % beta parameter of beta distribution
                  'tau_x_r', 0.002, ...               % rise time of feed-forward EPSP (s)
                  'tau_z_r', 0.002, ...               % rise time of lateral EPSP (s)
                  'tau_x_f', 0.02, ...                % fall time of feed-forward EPSP (s)
                  'tau_z_f', 0.02, ...                % fall time of lateral EPSP (s)
                  'tau_rf', 0.005, ...                % fall time of refractory EPSP (s)
                  'w_rf', -10, ...                    % refractory amplitude
                  'rec_delay', 0.005, ...             % delay on lateral synapses (s)
                  'eta', 0.005', ...                  % learning rate
                  'use_iw', false, ...                % use importance sampling
                  'save_interval', 1000, ...          % number of iterations between 2 save files
                  'reset_psps', false, ...            % reset EPSPs to 0 at beginning of sequence
                  'num_groups', 1, ...                % number of WTA circuits in population
                  'num_train_sets', 20000, ...        % number of training iterations
                  'pretrain_ff_weights', false, ...   % pretrain feed-forward weights
                  'use_variance_tracking', false, ... % use variance tracking
                  'update_on_spike', true, ...        % update right after spike event
                  'regen_patterns', true, ...         % regenerate patterns in each epoch
                  'regen_seqs', true, ...             % regenerate sequences in each epoch
                  'iw_track_speed', 0.001, ...        % speed of importance weight tracking
                  'num_samples', 1, ...               % number of samples at once
                };
            
    snn_options( 'colors', [ 0.0, 0.0, 1.0; ... %A
                 0.0, 0.8, 0.0; ... %B
                 1.0, 0.4, 0.4; ... %C
                 0.8, 0.8, 0.0; ... %D
                 0.7, 0.7, 0.7; ... %hold
                 0.8, 0.0, 0.8; ... %a
                 0.0, 0.6, 0.8; ... %b
                 0.8, 0.2, 0.2; ... %c
                 0.3, 0.6, 0.0; ... %d
                 0.4, 0.4, 0.4; ] );
             
    snn_options( 'FontName', 'Helvetica' );
    snn_options( 'FontSize', 8 );
end
