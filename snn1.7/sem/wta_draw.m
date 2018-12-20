function Z = wta_draw( P )
%
% Z = wta_draw( P )
%
% Draw a random sample vector Z from a distribution P.
% P is assumed to be normalised over columns, so for
% each column a sample is drawn. Output has the same
% dimensionality as P.
%
% 28.05.2010
%

    P_size = size(P);
    P_size(1) = 1;
    Z = cumsum( cumsum(P,1) >= repmat(rand(P_size),[size(P,1),ones(1,ndims(P)-1)]) ) == 1;
end
