function net = snn_alloc( net, allocators, reset )
% Allocate the fields defined by allocators.
%
% net = snn_alloc( net, allocators )
% net = snn_alloc( net, allocators, reset )
%
% Allocate the fields defined in the allocators cell-array.
%
% Generator function:
%   rand[A]     Random equally distributed numbers
%               between zero and A.
%   rand[A,B]   Random equally distributed numbers
%               between A and B.
%   randn[A,B]  Normal distributed numbers with
%               mean A and standard deviation B.
%   beta[A,B]   Beta distributed numbers.
%   spikes[A]   Boolean values, poisson distributed
%               with rate A.
%   zeros       Zero float values with double precision.
%   szeros      Zero float values with single precision.
%   izeros      Zero int32 values.
%   ones        Ones float values with double precision.
%   const[A]    Const float values with double precision
%               with value A.
%
% 17.11.2010
%

    if (nargin < 2)
        reset = false;
    end

    for i=1:3:length(allocators)
        
        field_name = allocators{i};
        
        gen_fct = allocators{i+1};
        
        [fct_args,a_pos] = snn_dispatch_args( net, gen_fct );
        
        gen_fct = gen_fct(1:a_pos-1);

        dim = snn_dispatch_args( net, allocators{i+2} );
        
        if isfield( net, field_name ) && ~reset
            if any( size(net.( field_name )) ~= dim )
                net.( field_name ) = alloc( fct_args, dim );
            end
        else
            net.( field_name ) = alloc( fct_args, dim );
        end
    end
    
%% local functions    
    function M = alloc( args, dim )
    % allocate an array of given size
        
        switch gen_fct
            case 'rand'
                if length(args)==1
                    M = args(1)*rand( dim );
                elseif length(args)==2
                    M = args(1) + (args(2)-args(1))*rand( dim );
                else
                    error('wrong number of arguments')
                end

            case 'randn'
                M = args(2)*randn( dim ) + args(1);
                
            case 'beta'
                M = betarnd( args(1), args(2), dim );
                
            case 'scbeta'
                M = arg(1)*betarnd( args(2), args(3), dim );
            
            case 'spikes'
                M = rand( dim )<=args(1);

            case 'szeros'
                M = zeros( dim, 'single' );
                
            case 'izeros'
                M = zeros( dim, 'int32' );
                
            case 'zeros'
                M = zeros( dim );
                
            case 'ones'
                M = ones( dim );

            case 'const'
                M = args(1)*ones( dim );
        end
    end
end
