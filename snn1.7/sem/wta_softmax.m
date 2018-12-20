function [P,A] = wta_softmax( U )
% 
% [P,A] = wta_softmax( U )
%
% Calculates the softmax function normalised
% over first dimension.
% Safe softmax function, avoids problems due to
% finite float precision.
%
% 02.06.2010
%


   U_max =  max(U(:),[],1);

   U_exp = exp( U - repmat( U_max, size(U) ) );
   
   A_n = sum(U_exp,1);
   
   P = U_exp./repmat( A_n, [size(U,1), ones(1,ndims(U)-1)] );
   
   if nargout>1
       A = log( A_n ) + U_max;
   end   
end
