function pcode = population_conding( Y, num_bins, min_max_range )
%
% Generate population coding of Y.
%

    if nargin<3
        min_max_range = minmax(Y);
    end

    [N,M] = size(Y);
    
    if size(num_bins,2)==1
        num_bins = repmat(num_bins,1,N);
    end
    
    if size(min_max_range,1)==1
        min_max_range = repmat(min_max_range,N,1);
    end
    
    num_bins = num_bins';
    
    pcode = sparse( sum(num_bins), M );

    for i=1:N
        
        Y_i = (Y(i,:)-min_max_range(i,1))/(min_max_range(i,2)-min_max_range(i,1));
    
        Y_i = max( Y_i, zeros(1,M) );
        Y_i = min( Y_i, ones(1,M) );
        Y_i = floor(Y_i*(num_bins(i)-1))+1;
    
        pcode(sum(num_bins(1:i-1))+1:sum(num_bins(1:i)),:) = ...
            sparse( Y_i, 1:M, ones( M, 1 ), num_bins(i), M );
    end
end
