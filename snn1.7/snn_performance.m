function performance = snn_performance( net, data )
% snn_performance: get performance measure for given net
%
% performance = snn_performance( net, data )
%
% Get performance measure for given net. The performance
% is calculated using the networks performance method.
%
% input
%   net:         A SNN network as created by
%                <a href="matlab:help snn_new">snn_new</a>.
%   data:        A test data structure or string pointing
%                to data file.
%
% output
%   performance: The performance measure calculated over the
%                whole data set.
%
% see also
%   <a href="matlab:help snn_new">snn_new</a>
%   <a href="matlab:help snn_list_methods">snn_list_methods</a>
%
% David Kappel 12.05.2010
%

    if ( nargin < 2 )        
        error( 'Not enough input arguments!' );
    end
    
    if ~isstruct( net )
        error( 'Unexpected argument type for ''net''!' );
    end

    if ischar( data )
        data = snn_load_data( data );
    end
    
    if ~isstruct( data )
        error( 'Unexpected argument type for ''data''!' );
    end
    
    performance = 0;
    num_perf = 0;
    
    for j=1:length(data)

        [net,perf_j] = net.p_performance_fcn( net, data(j) );
        
        if isfinite(perf_j)
            performance = performance + perf_j;
            num_perf = num_perf+1;
        end
    end

    if (num_perf>0)
        performance = performance/num_perf;
    else
        performance = nan;        
    end
end
