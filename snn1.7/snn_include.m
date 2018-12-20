function snn_include( varargin )
% snn_include: adds the given directories to the SNN search dirs
%
% snn_include
% snn_include( dir )
% snn_include( dir1, dir2, ... )
%
% Add the list of directories to the SNN search dirs. The
% search dirs are used to locate SNN methods. Calling
% snn_include without arguments will just add the SNN
% toolbox path to the matlab paths. The seach dirs can be
% shown by calling snn_options with argument 'search_dirs'
%
% see also
%   <a href="matlab:help snn_list_methods">snn_list_methods</a>
%   <a href="matlab:help snn_options">snn_options</a>
%
%
% David Kappel
% 17.11.2010
%
%

    if ( nargin < 1 )        
        addpath( pwd );
        return;
    end
    
    if ~iscellstr( varargin )
        
        if ischar( varargin )
            varargin = cellstr( varargin );
        else
            error( 'Unexpected argument type!' );
        end
    end

    search_dirs = snn_options( 'search_dirs' );
    
    if isempty( search_dirs )
        search_dirs = {};
    end

    for i = 1:length( varargin )
        
        dir = varargin{i};
    
        if ~ischar( dir )
            error( 'Unexpected argument type!' );
        end

        if ~exist( dir, 'dir' )
            
            home_path = what( [ 'snn', snn_version ] );
            
            if isempty( home_path )
                error( 'Toolbox home path not found!' );
            end
            
            dir = fullfile( home_path.path, dir );
            
            if ~exist( dir, 'dir' )            
                error( 'Directory ''%s'' not found!', dir );
            end
        end

        path_info = what( dir );

        if isempty( strmatch( path_info.path, search_dirs, 'exact' ) )

            search_dirs = { search_dirs{:} path_info.path };
            addpath( path_info.path );
        end
    end
    
    snn_options( 'search_dirs', search_dirs );
end
