function explore_data_set( data_set, i )
% Open figure to explore data set folder.
%
% 10.10.2011
% David Kappel
% 

    if (nargin < 2)
        i = 1;
    end
    
    fig = figure();
    
    index_edit = ...
        uicontrol( 'Style', 'edit', ...
                   'Parent', fig, ...
                   'String', num2str(i), ...
                   'Value', 0, ...
                   'UserData', data_set, ...
                   'Callback', @(hObject, eventdata)(index_edit_callback(hObject, eventdata)),...
                   'Position', [40 10 50 20]);
    
    uicontrol( 'Style', 'pushbutton', ...
               'Parent', fig, ...
               'String', '<',...
               'UserData', -1,...
               'Callback', @(hObject, eventdata)(next_button_callback(hObject, eventdata)),...
               'Position', [20 10 20 20] );

    uicontrol( 'Style', 'pushbutton', ...
               'Parent', fig, ...
               'String', '>',...
               'UserData', 1,...
               'Callback', @(hObject, eventdata)(next_button_callback(hObject, eventdata)),...
               'Position', [90 10 20 20] );

    index_edit_callback( index_edit, [] );
           
    % Callback function for 'index' edit field.
    function index_edit_callback(hObject, eventdata)

        cur_val = get( hObject, 'Value' );
        input = get( hObject, 'String' );
        
        new_val = str2double(input);
        
        if isnan( new_val )
            set( hObject, 'String', num2str(cur_val) );
        elseif ( cur_val ~= new_val )
            base_path = get( hObject, 'UserData' );
            subdirs = dir( [ base_path, 'data_set_*.mat' ] );
            
            if (new_val<1) || (new_val>length(subdirs))
                return
            end
            
            parent = get( hObject, 'Parent' );
            siblings = get( parent, 'Children' );
            set( hObject, 'Enable', 'off' );
            set( siblings(end-1), 'Enable', 'off' );
            set( siblings(end-2), 'Enable', 'off' );
            drawnow;            

            try                
                plot_data_set( base_path, subdirs(new_val).name, parent );
                set( hObject, 'Value', new_val );
                set( hObject, 'String', num2str(new_val) );
            catch
                set( siblings(end-3), 'String', 'Error while loading data set!' );
                set( hObject, 'String', num2str(cur_val) );
                set( hObject, 'Value', cur_val );
                
                fprintf( 'Error while loading data set!\n%s\n', lasterr );
            end
            
            set( hObject, 'Enable', 'on' );
            set( siblings(end-1), 'Enable', 'on' );
            set( siblings(end-2), 'Enable', 'on' );
            drawnow;
        end
    end
    
    % Callback function for 'next' button.
    function next_button_callback(hObject, eventdata)
        
        parent = get( hObject, 'Parent' );
        increment = get( hObject, 'UserData' );
        siblings = get( parent, 'Children' );
        cur_val = get( siblings(end), 'Value' );        
        set( siblings(end), 'String', num2str( cur_val+increment ) );
        index_edit_callback( siblings(end), [] );
    end

end
