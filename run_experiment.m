% run_experiment.m
% 
% Train stochastic trajectories into a recurrent
% network of SEMs.
%
% This setting realizes experiment 1 of:
%   D. Kappel, B. Nessler, and W. Maass. STDP Installs in
%   Winner-Take-All Circuits an Online Approximation to
%   Hidden Markov Model Learning. PLOS Computational
%   Biology, 2014.
%
%
% Institute for Theoretical Computer Science
% Graz University of Technology
%
% 10.10.2011
% David Kappel
% http://www.igi.tugraz.at/kappel/
%


%% init

run snn1.7/snn_include;
snn_include( 'sem', 'hmmsem', 'plotting' );

default_params = default_options();


%% train the network

ex_path = do_learning_task( 'results/', ...
                  { [ 1,2,5,6,7 ], ...     %AB-hold-ab
                    [ 2,1,5,7,6 ], ...     %BA-hold-ba
                    [ 3,4,5,8,9 ], ...     %CD-hold-cd
                    [ 4,3,5,9,8 ] }, ...   %DC-hold-dc
                  default_params{:}, ...
                  'free_run_seqs', { [ 1,2,5 ], ...    %AB-hold
                                     [ 2,1,5 ], ...    %BA-hold
                                     [ 3,4,5 ], ...    %CD-hold
                                     [ 4,3,5 ] }, ...  %DC-hold
                  'pat_labels', { 'A', 'B', 'C', 'D', 'hold', 'a', 'b', 'c', 'd' }, ...
                  'free_run_pat_lenghts', [0.050,0.050,0.150], ...
                  'num_neurons', 100, ...             % number of WTA neurons
                  'num_inputs', 200, ...              % number of afferent neurons
                  'free_run_time', 0.400, ...         % free run time (s)
                  'save_interval', 100, ...           % number of iterations between 2 save files
                  'num_train_sets', 5000, ...        % number of training iterations
                  'collect', '[At,R]', ...
                  'num_epochs', 1, ...
                  'changelog_flag', 'N' );


%% evaluate training result

fprintf('\n\nevaluating training perforamnce:\n')
              
eval_mem_task( ex_path );



%% display results

data = load( [ex_path,'mem_task_test.mat'] );

plot_mem_task( data, [], 'seq_id', 1, 'data_set', 'test_data', ...
               'plot_spikes', false, 'neuron_order', data.I, 'fig_file', 'evoced_1', ...
               'peth_set', 'peth_test', 'dy', 0.5, 'base_path', ex_path );

plot_mem_task( data, [], 'seq_id', 2, 'data_set', 'test_data', ...
               'plot_spikes', false, 'neuron_order', data.I, 'fig_file', 'evoced_2', ...
               'peth_set', 'peth_test', 'dy', 0.5, 'base_path', ex_path );

explore_data_set( ex_path )

