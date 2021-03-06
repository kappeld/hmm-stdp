function [Z,k] = wta_draw_k( P )
%
% [Z,k] = wta_draw_k( P )
%
% Draw a random sample vector Z from a distribution P.
% P is assumed to be normalised over columns, so for
% each column a sample is drawn. Output has the same
% dimensionality as P.
%
% 05.08.2011
%

    k = find( cumsum(P,1) > rand(), 1 );
    
    if isempty(k)
        k = size(P,1);
    end
    
    Z = sparse( k, 1, 1, size(P,1), 1 );    
end
