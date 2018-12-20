function [net, performance] = snn_performance_ct( net, data )
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

    testTprob = double( min( 1, data.Pt(1:(end-1),:) ) );
    
    if isfield(data, 'T')
        T = sparse( double( data.T ), 1:length(double( data.T )), 1 );
        targets = T(:,ceil(data.Pt(end,:)*data.sample_rate));
    else
        error('FIXME: Not implemented yet!!!');
    end
    
    probY = sum(targets,1);
    probLT = sparse(targets)*diag(probY)*testTprob';
    probLT = probLT/sum(sum(probLT));
    
    HLT=sum(sum(probLT.*log(max(eps,probLT))));
    HT=sum(sum(probLT,1).*log(max(eps,sum(probLT,1))));
    HL=sum(sum(probLT,2).*log(max(eps,sum(probLT,2))));

    performance = (HLT-HT)/HL;
end
