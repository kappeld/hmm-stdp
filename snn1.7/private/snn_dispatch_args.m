function [args,start_pos,end_pos] = snn_dispatch_args( net, token )
% Dispatch an argument string
%
% [args,start_pos,end_pos] = snn_dispatch_args( net, token )
%
% Dispatch an rgument of the form [<arg1>,<arg2>,...,<argN>].
% The argument fields can be any real valued argument, or a
% field name referencing to a field inside the net structure.
% snn_dispatch_args always returns a matrix containing real
% values or generates an error if any of the arguments inside
% the token string is eighter not real valued or not referencing
% to a field.
%
% 17.11.2010
%

    [toks,start_pos,end_pos] = snn_parse_args( token );

    args = zeros( 1, numel(toks) );

    for d = 1:length( toks )

        value = str2double( toks{d} );

        if strcmp( toks{d}, 'nan' )
            args(d) = nan;
        elseif ~isnan( value )
            args(d) = value;
        else

            if ~isfield( net, toks{d} )
                error( 'Parameter ''%s'' undefined in allocator comamnd ''%s''!', ...
                       toks{d}, token );
            end
            
            args(d) = net.(toks{d});
        end
    end
end
