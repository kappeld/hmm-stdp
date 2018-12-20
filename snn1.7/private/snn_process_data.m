function [net, data_out] = snn_process_data( net, fcn_data_processor, data_gen, ...
                                             step_size, num_runs, keep_data_order, ...
                                             lead_in, field_collectors )
% snn_process_data: process the data set using the given processor function
%
% [net, data_out] = snn_process_data( net, fcn_data_processor, data, ...
%                                     step_size, num_runs, keep_data_order, ...
%                                     lead_in, field_collectors )
%
% Process the data set using the given processor function. The
% function handler fcn_data_processor is called for each data
% block of given data set.
%
% input
%   net:                A SNN network.
%                       See <a href="matlab:help snn_new">snn_new</a>.
%   fcn_data_processor: Function handler pointing to a function with interface:
%                       [net,Z,P] = fcn_data_processor( net, data, ct )
%   data:               A snn data structure or a string containing a file
%                       patterns that points to files from which the data
%                       structures should be loaded.
%                       See <a href="matlab:help snn_load_data">snn_load_data</a>.
%   step_size:          The block size of the data blocks.
%   num_runs:           The number of times the data set should be processed.
%   keep_data_order:    If true the data is processed in the given order, if
%                       false the order of the data is randomised before
%                       being processed.
%   lead_in:            An integer number. snn_process_data patches the data
%                       with the given number of time steps from previous
%                       data set.
%   field_collectors:   A cell array strings that holds the field names of
%                       fields that are collected from the net structure.
%                       The fields are collected after each block that was
%                       processed and stored in the data_out structure.
%
% output
%   net:                The modified network structure.
%   data_out:           The processed data containing the default fields X
%                       and P and any additional fields that where
%                       collected by the field_collectors.
%
% David Kappel 6.12.2010
%

    if keep_data_order
        perm = 1:num_runs;
    else
        perm = randperm(num_runs);
    end
    
    for i=1:num_runs

        data_set = data_gen.fcn_generate( data_gen, perm(i) );
        
        if isfield( net, 'reset_psps' ) && net.reset_psps
            if isfield( net, 'hX_init' )
                net.hX = net.hX_init;
                net.hZ = net.hZ_init;
            else
                net.hX(:) = 0;
                net.hZ(:) = 0;
            end
            
            if isfield( net, 'last_spike_t' )
                net.last_spike_t(:) = 0;
            end
        end

        [ net, ...
          data_set.Zt, ...
          data_set.Pt ] = fcn_data_processor( net, data_set, [1,length(data_set.time)] );

        for c = 1:length(field_collectors)
            data_set.(field_collectors{c}) = net.(field_collectors{c});
        end

        if ( nargout > 1 )
            data_out(i) = data_set;
        end

        if isfield(net, 'p_user_interupt')
            break;
        end
    end
end
