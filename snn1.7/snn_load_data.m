function data = snn_load_data( file_pattern, varargin )
% snn_load_data: loads experiment data form disk
%
% data = load_experiment( file_pattern )
% data = load_experiment( file_pattern, ... )
%
% Loads experiment data form disk. Additional to the file_pattern
% an arbitrary number of field names can be given that should be
% loaded. If no field names are given all fields are read form
% the fieles.
%  
% A SNN data file is a .mat file that (at least) contains the fields:
%     X: network input.
%
% input
%   file_pattern:  A string containing a file name or a pattern containing
%                  placeholders like '*' if multiple files should be loaded.
%
% output
%   data:          The SNN data structure containing the loaded fields.
%
% David Kappel 07.06.2010
%

    if (nargin<1)
        errro( 'Not enough input arguments!' );
    end

    files = dir( file_pattern );
    filedir = fileparts( file_pattern );
    
    verbose = snn_options( 'verbose' );
    
    if isempty(verbose)
        verbose = true;
    end

    if isempty(files)
        error( 'No files found that match pattern ''%s''!', file_pattern );
    end
    
    for i=1:length(files)
        if verbose
            disp( [ '  loading file : ', files(i).name ] );
        end
        data(i) = load( fullfile( filedir, files(i).name ), varargin{:} );
    end
end
