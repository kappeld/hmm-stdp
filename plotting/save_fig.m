function save_fig( h_ref, file_name, base_path )
% Saves the figure given by h_ref to the specified file_name.
% save_fig creates 3 files in the 'figures' subfolder storing
% the figure in fig, eps and pdf format with 300 dpi
% resolution.
%
% 30.09.2011
% David Kappel
% 

    if nargin < 3
        base_path = '';
    end

    file_path = [ base_path , 'figures/', file_name ];
    
    mkdir([ base_path , 'figures/']);

    set( h_ref, 'PaperPositionMode', 'auto' );
    
    f_h = [ '-f', num2str(h_ref) ];

    saveas( h_ref, file_path, 'fig' );
    print( f_h, '-depsc', '-painters', '-r300', file_path );
    print( f_h, '-dpdf', '-painters', '-r300', file_path );
    
    fprintf('figure saved to ''%s''\n',file_path)
end
