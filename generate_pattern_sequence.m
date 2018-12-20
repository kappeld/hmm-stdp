function data = generate_pattern_sequence( patterns, pat_sequences, ...
                                           pattern_length, num_sequences, ...
                                           targets_per_pattern, length_std, pat_labels )
% generate_pattern_sequence
%
% Generates a spike train sequence composed of rate patterns of a given length.
%
% 05.11.2011
% David Kappel
%

    seq_id = [];
    
    time_padding = 0;
    
    st_process = 'poisson';
                                       
    if (nargin <= 2) && isstruct( patterns )
        set_generator = patterns;
        patterns = set_generator.patterns;
        seq_id = set_generator.seq_id;
        pattern_length = set_generator.pattern_length;
        targets_per_pattern = set_generator.targets_per_pattern;
        length_std = set_generator.length_std;
        num_sequences = 1;
        
        if isfield( set_generator, 'process' )
            st_process = set_generator.process;
        end
        
        if isfield( set_generator, 'time_padding' )
            time_padding = set_generator.time_padding;
        end
        
        if (nargin == 2) && ischar( pat_sequences ) && strcmp( pat_sequences, 'free_run' )
            if isfield( set_generator, 'free_run_seqs' )
                pat_sequences = set_generator.free_run_seqs{ seq_id };
            else            
                pat_sequences = set_generator.pat_sequences{ seq_id };
                pat_sequences(end) = length( patterns );
            end
        else
            pat_sequences = set_generator.pat_sequences{ seq_id };
        end
        
    elseif (nargin < 5)
        error( 'Not enought input arguments!' );
    elseif (nargin < 6)
        length_std = 0;
    end
    
    num_inputs = size( patterns{1}, 1 );
    
    num_pats = length( patterns );
    
    if iscell( pat_sequences )        
        num_patterns = length( pat_sequences );
    else
        [ num_patterns, pats_per_input ] = size( pat_sequences );
     	pat_sequences = mat2cell( pat_sequences, 1, repmat( pats_per_input, 1, num_patterns ) );
    end
    
    if ( numel( pattern_length ) == 1 )
        pattern_length = repmat( pattern_length, length( pat_sequences{1} ), 1 );
    end
    
    if ( numel( length_std ) == 1 )
        length_std = repmat( length_std, length( pat_sequences{1} ), 1 );
    end
        
    if (nargin < 7)
        pat_labels = mat2cell( char( 'A' + (1:num_pats) - 1 ), 1, ones(1,num_pats) );
    end
    
    sample_rate = 1000;
    
    data = struct();
    
    for i = 1:(num_sequences)

        seq = pat_sequences{ mod(i-1,num_patterns)+1 };
        pats_per_input = length( seq );
        
        pat_lengths = zeros(pats_per_input,1);

        Xt = [];
        Lt = [];
        time = 0;
        labels = struct;    

        t = 1;
        for l=1:pats_per_input

            pat_no = seq(l);

            pat_lengths(l) = randn()*length_std(l) + pattern_length(l);
            
            if (pat_lengths(l) < 0.010)
                pat_lengths(l) = 0.010;
                warning( 'Minimum pattern length is 10ms' );
            end
            
            pat = generate_rate_pattern( patterns{ pat_no }, pat_lengths(l), st_process );
            
            pat(2,:) = pat(2,:) + sum( pat_lengths(1:(l-1)) );
            
            num_spikes = size(pat,2);

            Xt = [ Xt, pat ];
            Lt = [ Lt, [pat(1,:); repmat(pat_no, 1, num_spikes); pat(2,:)] ];
            time = [ time, pat(2,:) ];

            labels(l).descriptor = pat_labels{pat_no};
            labels(l).start_sample = t;
            labels(l).stop_sample = (t+num_spikes);

            t = t+num_spikes;
        end
        
        time = [ time, sum( pat_lengths )+time_padding ];

        total_length = ceil( time(end)*sample_rate );

        T = ceil( targets_per_pattern*pats_per_input*((1:total_length)/total_length));

        data(i).Xt = Xt;
        data(i).Lt = Lt;
        data(i).T = T;

        data(i).time = time;
        data(i).x_range = single( 1:num_inputs );
        data(i).labels = labels;
        data(i).sample_rate = sample_rate;
        data(i).seq_id = seq_id;
    end
end
