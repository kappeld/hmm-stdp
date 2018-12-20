function [varargout] = snn_options( varargin )
% snn_options: sets and gets global SNN options
%
% options = snn_options
% value = snn_options( op_name )
% value = snn_options( op_name_1, op_value_1, op_name_2, op_value_2, ... )
%
% Sets and gets global SNN options. If op_value is given the option is set
% to the new value and the old one is returned. If the option does not
% exist, or is set the first time [] is returned. If no option name is
% provided the whole options structure is rectured.
%
% inputs
%   op_name:   A string containing the option name.
%
% David Kappel
% 20.08.2010
%

    global snn_options__;
    
    if isempty( snn_options__ )        
        snn_options__ = struct;        
    end

    if ( nargin > 0 )
        
        for k=1:2:length(varargin)
        
            op_name = varargin{k};
                        
            if ~ischar( op_name )
                error( 'Unexpected argument type for ''op_name''!' );
            end

            if isfield( snn_options__, op_name )
                varargout((k+1)/2) = { snn_options__.(op_name) };
            else
                varargout((k+1)/2) = { [] };
            end

            if ( nargin > k )
                
                new_value = varargin{k+1};
                
                if isempty( new_value )
                    if isfield( snn_options__, op_name )
                        snn_options__ = rmfield( snn_options__, op_name );
                    end
                else
                    snn_options__.(op_name) = new_value;
                end
            end
        end
    else        
        varargout(1) = { snn_options__ };
    end
end
