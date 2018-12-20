function [net,Z,P] = snn_sample_sid( net, data, ct )
% Draws a spike from sem network - self inhibition and recurrent delays
%
% [net,Z,P] = snn_sample_ct( net, data, ct )
%
% Generates the output of a given wta-network by drawing
% a winner neuron from the current output probabilities
% (standard version).
%
% inputs:
%   net:  A wta-network, see: wta-new()
%   data: A data structure to be simulated.
%   ct:   Current time index.
%
% output:
%   net:      The (modified) network structure
%   Z:        Network output spikes
%
% @parameters:
%   tau_x_r         0.002        time constant of epsp window rise
%   tau_z_r         0.002        time constant of epsp window rise
%   tau_x_f         0.02         time constant of epsp window fall
%   tau_z_f         0.02         time constant of epsp window fall
%   tau_rf          0.005        time constant of refractory window
%   w_rf            -10          amplitude of refractory window
%   lambda          2000         network spike rate   
%   mean_rec_delay  0.010        mean delay on recurrent synapses
%   std_rec_delay   0.005        standard deviation of recurrent delay
%   mean_ff_delay  0.010         mean delay on recurrent synapses
%   std_ff_delay   0.005         standard deviation of recurrent delay
%
% @fields:
%   rec_delay     randn[mean_rec_delay,std_rec_delay]  [num_neurons,num_neurons]
%   ff_delay      randn[mean_ff_delay,std_ff_delay]    [num_neurons,num_input]
%   rec_spikes    zeros [0,0]            recurrent spike events
%   last_spike_t  zeros [num_neurons,1]  last output spike times
%   At    szeros  [2,0]             log activity of network at spike times
%   hX    zeros   [num_inputs,2]    feedforward PSPs
%   hZ    zeros   [num_neurons,2]   recurrent PSPs
%   R     zeros   [1,1]             current importance weight
%   num_o zeros   [num_neurons,1]   number of output spikes per neuron
%   

    t = data.time(ct(1));
    
    time_range = data.time(ct(end))-data.time(ct(1));    
    num_spikes = poissrnd( double( net.lambda*time_range ) );
    spike_times = sort( time_range*rand( 1, num_spikes ) );    
    rec_spikes = [ net.rec_spikes, inf(2,num_spikes) ];
    num_rs = size(net.rec_spikes,2);
    
    Z = zeros(2, num_spikes, 'single');
    P = zeros(net.num_neurons+1, num_spikes, 'single');
    net.At = zeros(2,num_spikes, 'single');
    
    last_input_spikes = repmat(t, net.num_inputs, 1);
    last_output_spikes = repmat(t, net.num_neurons, 1);
    last_spike_times = net.last_spike_t;
    
    hX = net.hX;
    hZ = net.hZ;    
    eta = net.eta;
    
    W_exp = exp( double( net.W ) );
    V_exp = exp( double( net.V ) );
    V0_exp = exp( double( net.V0 ) );
    
    net.d_W =  zeros(net.num_neurons,net.num_inputs);
    net.d_V =  zeros(net.num_neurons,net.num_neurons);
    net.d_V0 = zeros(net.num_neurons,1);
    
    i = ones(net.num_neurons,1);
    l = ones(net.num_neurons,1);
    
    hX_all = zeros(net.num_inputs,num_spikes);
    hZ_all = zeros(net.num_neurons,num_spikes);
    
    try
        
        for j = 1:num_spikes

            t = spike_times(j);

            for k = 1:net.num_neurons
                while (i < size(data.Xt,2)) && (t > data.Xt(2,i))

                    n_id = data.Xt(1,i);
                    sp_t = data.Xt(2,i);

                    hX(n_id,1) = hX(n_id,1)*exp(-double(sp_t-last_input_spikes(n_id))/net.tau_x_r) + net.W(k,n_id);
                    hX(n_id,2) = hX(n_id,2)*exp(-double(sp_t-last_input_spikes(n_id))/net.tau_x_f) + net.W(k,n_id);

                    last_input_spikes(n_id) = sp_t;
                    i = i+1;
                end
            
                while (l < size(rec_spikes,2)) && (t > rec_spikes(2,l))

                    n_id = rec_spikes(1,l);
                    sp_t = rec_spikes(2,l);

                    hZ(n_id,1) = hZ(n_id,1)*exp(-double(sp_t-last_output_spikes(n_id))/net.tau_z_r) + net.V(k,n_id);
                    hZ(n_id,2) = hZ(n_id,2)*exp(-double(sp_t-last_output_spikes(n_id))/net.tau_z_f) + net.V(k,n_id);

                    last_output_spikes(n_id) = sp_t;
                    l = l+1;
                end
            end

            hZ(:,1) = hZ(:,1).*exp(-double(t-last_output_spikes)/net.tau_z_r);
            hZ(:,2) = hZ(:,2).*exp(-double(t-last_output_spikes)/net.tau_z_f);
            hX(:,1) = hX(:,1).*exp(-double(t-last_input_spikes)/net.tau_x_r);
            hX(:,2) = hX(:,2).*exp(-double(t-last_input_spikes)/net.tau_x_f);

            d_hX = diff(hX,1,2);
            d_hZ = diff(hZ,1,2);
            
            hX_all(:,j) = d_hX;
            hZ_all(:,j) = d_hZ;
            
            U = d_hX + d_hZ + net.V0 - exp(-double(t-last_spike_times)/net.tau_rf);
            [P_t,A] = wta_softmax( U );

            [Z_i,k] = wta_draw_k(P_t);

            Z(:,j) = [k;t];
            P(:,j) = [P_t;t];
            
            last_input_spikes(:) = t;
            last_output_spikes(:) = t;
            last_spike_times(k) = t;
            
            rec_spikes(:,j+num_rs) = [k,t+net.rec_delay(k)];
            
            net.num_o(k) = net.num_o(k)+1;

            net.d_W(k,:) = net.d_W(k,:) + (d_hX'-W_exp(k,:))./max(eta,W_exp(k,:));
            net.d_V(k,:) = net.d_V(k,:) + (d_hZ'-V_exp(k,:))./max(eta,V_exp(k,:));
            net.d_V0     = net.d_V0     + (Z_i-V0_exp)./max(eta,V0_exp);

            net.At(:,j) = [A,t];
        end
    
    catch
        fprintf('There has been an error, while sampling!\nExcluding run from training\n');
        file_name = sprintf('/tmp/error_data_%u_%04u.mat', net.iteration, round(9999*rand()) );
        the_error = lasterror();
        fprintf( '\n%s\n', the_error.message );
        fprintf( '  in %s on line %i\n\n', the_error.stack(1).name, the_error.stack(end).line );
        fprintf('saving workspace to: %s\n', file_name);
        save( file_name );
        net.R = -100000;
        fprintf('pausing 10 seconds...\n');
        pause(10);
        return;
    end
    
    
    net.rec_spikes = rec_spikes(:,l:end);
    net.rec_spikes(2,:) = net.rec_spikes(2,:) - t;
    net.last_spike_t = last_spike_times - t;
    
    net.R = sum( net.At(1,:) );    
    net.hX = hX;
    net.hZ = hZ;
    net.hX_all = hX_all;
    net.hZ_all = hZ_all;

end
