function snn_assert_equal( val1, val2, tolerance )
% snn_assert_equal: assert the two values to be equal
%
% snn_assert_equal( val1, val2 )
% snn_assert_equal( val1, val2, tolerance )
%
% Compares the two values and raises an error message
% if their eucledian distance is above the given
% tolerance. If no tolerance is given, the value default
% value of 1e-10 is taken.
%
% see also
%   <a href="matlab:help snn_options">snn_options</a>
%
%
% David Kappel
% 12.12.2010
%
%

    if ( nargin < 2 )
        error( 'Not enought input arguments!' );
    end
    
    if ( nargin < 3 )
        tolerance = 1e-10;
    end
    
    if sum( ( val1 - val2 ).^2 ) > tolerance
        error( 'Assertion failed!' );
    end
end
