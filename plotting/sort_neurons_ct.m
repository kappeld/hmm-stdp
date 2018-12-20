function idx = sort_neurons_ct( data, field_name )
%
% Generate the neural sorting index.
%

    P = [];
    
    for i=1:length(data)
        P = [ P, data(i).Pt(1:end-1,:) ];
    end

    [v,a] = max(P,[],2);
    [v,idx] = sort(a);
end
