function [net, performance] = snn_performance_ce( net, data )
% Normalized conditional entropy performance measure.
%
% [net, performance] = wta_performance_ce( net, data )
%
% SNN performance method that calculates the normalized
% conditional entropy.
%
% net:
%   A wta-network, see wta_new().
%
% data:
%   Test data structure.
%

    testTprob = min( 1, data.P );

    probY = sum(data.T,1);
    probLT = sparse(data.T)*diag(probY)*testTprob';
    probLT = probLT/sum(sum(probLT));
    
    HLT=sum(sum(probLT.*log(max(eps,probLT))));
    HT=sum(sum(probLT,1).*log(max(eps,sum(probLT,1))));
    HL=sum(sum(probLT,2).*log(max(eps,sum(probLT,2))));

    performance = (HLT-HT)/HL;
end
