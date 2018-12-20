function file_name = locate_data_file( data_set_dir, i, epoch )
% Locate the data set file in the given directory
%
% file_name = locate_data_file( data_set_dir, index )
%
% Locate the data set file with given index in the given
% directory and returns the file location. data_set_dir
% must be a string pointing to a directory, index must
% be an interger. If no matching file was found, an empty
% string is returned.
%

    if (nargin < 3)
        epoch = [];
    end

    if  ~isempty(data_set_dir) && ( data_set_dir(end) ~= '/' )
        data_set_dir(end+1) = '/';
    end
    
    if ischar( i )
        if exist( [ data_set_dir, i ], 'file' )
            file_name = [ data_set_dir, i ];
        elseif strcmp(i,'end') || strcmp(i,'last')
            file_names_dir = dir( [ data_set_dir, 'data_set_*.mat' ] );
            file_name = [ data_set_dir, file_names_dir(end).name ];
        elseif strcmp(i,'first')
            file_names_dir = dir( [ data_set_dir, 'data_set_*.mat' ] );
            file_name = [ data_set_dir, file_names_dir(1).name ];
        else
            error( 'file not found ''%s''', [ data_set_dir, i ] );
        end
    else
        if isempty(epoch)
            file_names_dir = dir( [ data_set_dir, sprintf( 'data_set_*%d.mat', i ) ] );
        else
            file_names_dir = dir( [ data_set_dir, sprintf( 'data_set_*%d%05d.mat', epoch-1, i ) ] );
        end

        file_name_len = inf;

        file_name = '';

        for i = 1:length( file_names_dir )
            if  ( length( file_names_dir(i).name ) < file_name_len )
                file_name = [ data_set_dir, file_names_dir(i).name ];
                file_name_len = length( file_names_dir(i).name );
            end
        end
    end
end
