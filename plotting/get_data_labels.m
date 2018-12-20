function labels = get_data_labels( Zt, pat_labels, init_lbls )
%
% Get the data lables from the porvided network output.
%

    window_length = 0.010;
    min_win_length = 0.020;

    if (nargin > 2)
        labels = init_lbls;
        if isempty( fieldnames( init_lbls ) )
            l_offs = 0;
        else
            l_offs = length(init_lbls);
        end
    else    
        labels = struct();
        l_offs = 0;
    end
    
    num_patterns = max( Zt(2,:) );
    
    vote_map = zeros( num_patterns, 1 );
    
    win_t = Zt(3,1);
    start_sample = 1;
    
    last_spike = 1;
    
    current_vote = -1;
    
    i = 1;
    
    for t = 1:size(Zt,2)

        cur_t = Zt(3,t);
        
        if ( Zt(2,t) > 0 )
            vote_map( Zt(2,t) ) = vote_map( Zt(2,t) ) + 1;
        end
        
        while ( cur_t - win_t > window_length )
            if ( Zt(2,last_spike) > 0 )
                vote_map( Zt(2,last_spike) ) = vote_map( Zt(2,last_spike) )-1;
            end
            last_spike = last_spike + 1;            
            win_t = Zt(3,last_spike);
        end
        
        [support,vote] = max( mk_normalised( vote_map ) );
        
        if ( support < 0.5 )
            vote = -1;
        end
        
        if ( current_vote ~= vote )
            
            win_length = Zt(3,t) - Zt(3,start_sample);
            
            if ( win_length > min_win_length  ) && ( current_vote > 0 )
                
                labels(i+l_offs).start_time = max(0,Zt(3,start_sample)-window_length/2);
                labels(i+l_offs).stop_time = max(0,Zt(3,t)-window_length/2);
                labels(i+l_offs).descriptor = pat_labels{current_vote};
                labels(i+l_offs).show_border = true;
                i = i+1;
            end
            
            current_vote = vote;
            
            start_sample = t;            
        end        
    end
    
    if (i==1)
        labels = [];
    end
    
    win_length = Zt(3,t) - Zt(3,start_sample);

    if ( win_length > min_win_length  ) && ( current_vote > 0 )

        labels(i).start_time = Zt(3,start_sample);
        labels(i).stop_time = Zt(3,t);
        labels(i).descriptor = pat_labels{current_vote};
    end
end
