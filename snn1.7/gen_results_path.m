function save_path = changelog_gen_path( exp_path, flag )
% generate path within the given directory
%
% save_path = changelog_gen_path( path )
%
% David Kappel
% 10.10.2011
%

    cur_time = round(clock);

    exp_path_dirs = dir( exp_path );
    
    cur_save_path_idx = 0;

    for f_i = length( exp_path_dirs ):-1:1
        if exp_path_dirs( f_i ).isdir
            cur_save_path_idx = f_i;
            break;
        end
    end

    if cur_save_path_idx
        
        cur_save_path = exp_path_dirs( cur_save_path_idx ).name;
        
        fprintf( 'Last session found: %s.\n', cur_save_path );
        if (nargin<2) || isempty(flag)
            user_input = input( 'Start new session (N), Continue last session (c), Continue other session (o), Testing (t), Abord (a): ', 's' );
        else
            user_input = flag;
        end
        
        if lower( user_input ) == 'a'
            save_path = [];
            return;
        elseif lower( user_input ) == 'c'
            save_path = [ exp_path, cur_save_path, '/' ];
        elseif lower( user_input ) == 'o'
            
            path_idxs = nan(1,length( exp_path_dirs ));
            num_paths = 0;
            
            for f_i = length( exp_path_dirs ):-1:1
                if exp_path_dirs( f_i ).isdir && ( exp_path_dirs( f_i ).name(1) ~= '.' )
                    num_paths = num_paths+1;
                    path_idxs(num_paths) = f_i;
                end
            end            
            
            fprintf( 'Sessions found:\n' );
            
            for f_i=1:num_paths
                fprintf( '  [%2i] %s', f_i-1, exp_path_dirs( path_idxs(f_i) ).name );
                if (mod(f_i,4) == 0)
                    fprintf( '\n' )
                end
            end
            
            u_input = input( '\nSelect session: ', 's' );
            
            if ~isempty( str2num( u_input ) )
                cur_save_path = exp_path_dirs( path_idxs(str2num( u_input )+1) ).name;
            end
            
            save_path = [ exp_path, cur_save_path, '/' ];

        elseif lower( user_input ) == 't'
            save_path = [ tempname, '/' ];
            mkdir( save_path );
            fprintf( 'Saving session to temporary path: %s.\n', save_path );
            return;
        else
            save_path = [ exp_path, sprintf( '%04i.%02i.%02i-%02ih%02im%02is/', cur_time ) ];
        end
        
    else
        save_path = [ exp_path, sprintf( '%04i.%02i.%02i-%02i:%02i:%02i/', cur_time ) ];
    end
    
    fprintf( 'Saving session data to: %s.\n', save_path );
    
    changelog_file = 'changelog.mat';
    source_path = [ save_path, 'source/' ];
    
    mkdir( save_path );

end
