function pat = generate_rate_pattern( pattern, pattern_length, process )
% generate_rate_pattern
%
% Generates a rate patterns of a given length.
%
% 05.11.2011
% David Kappel
%

    if (nargin < 3)
        process = 'poisson';
    end
    
    switch process
        case 'poisson'
            num_spikes = poissrnd( pattern*pattern_length );

            pat = zeros(2,sum(num_spikes),'single');

            start_i = cumsum( [1; num_spikes(1:end-1)] );

            for i = 1:length(num_spikes);

                pat( 1, start_i(i):(start_i(i)+num_spikes(i)-1) ) = i;
                pat( 2, start_i(i):(start_i(i)+num_spikes(i)-1) ) = ...
                    pattern_length*rand( 1, num_spikes(i) );
            end

            [v,idx] = sort( pat(2,:) );
            pat = pat(:,idx);
            
        case 'regular'
            num_spikes = round( pattern*pattern_length );

            pat = zeros(2,sum(num_spikes),'single');

            start_i = cumsum( [1; num_spikes(1:end-1)] );

            for i = 1:length(num_spikes);

                pat( 1, start_i(i):(start_i(i)+num_spikes(i)-1) ) = i;
                pat( 2, start_i(i):(start_i(i)+num_spikes(i)-1) ) = pattern_length*(0:(1/num_spikes(i)):(1-(1/num_spikes(i))));
            end

            [v,idx] = sort( pat(2,:) );
            pat = pat(:,idx);
    end
end
