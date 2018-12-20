function [net,Z,P] = snn_sample_hmmrf( net, data, ct )
% Draws a spike from sem network - continuous time spiking hmm sems realistic refractory mechanim.
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
%   w_rf            -10          strength of refractroy window
%   lambda          2000         network spike rate
%   mean_rec_delay  0.010        mean delay on recurrent synapses
%   std_rec_delay   0.000        standard deviation of recurrent delay
%   use_variance_tracking  false
%
% @fields:
%   rec_delay     randn[mean_rec_delay,std_rec_delay]  [num_neurons,1]
%   rec_spikes    zeros [0,0]              recurrent spike events
%   last_spike_t  zeros  [num_neurons,1]   last output spike times
%   At    szeros  [2,0]             log activity of network at spike times
%   hX    zeros   [num_inputs,2]    feedforward PSPs
%   hZ    zeros   [num_neurons,2]   recurrent PSPs
%   R     zeros   [1,1]             current importance weight
%   num_o zeros   [num_neurons,1]   number of output spikes per neuron
%   eta_W const[eta] [num_neurons,num_inputs]
%   eta_V const[eta] [num_neurons,num_neurons]
%   eta_0 const[eta] [num_neurons,1]
%   SW    zeros      [num_neurons,num_inputs]
%   QW    ones       [num_neurons,num_inputs]
%   SV    zeros      [num_neurons,num_neurons]
%   QV    ones       [num_neurons,num_neurons]
%   S0    zeros      [num_neurons,1]
%   Q0    ones       [num_neurons,1]
%   

    t = data.time(ct(1));
    
    time_range = data.time(ct(end))-data.time(ct(1));
    
    %num_spikes = round(net.lambda*time_range)-1;
    %d_spikes = (time_range/(num_spikes));
    %spike_times = data.time(ct(1))+d_spikes:d_spikes:data.time(ct(end));
    
    %num_spikes = round(net.lambda*time_range);
    num_spikes = poissrnd( double( net.lambda*time_range ) );
    spike_times = sort( time_range*rand( 1, num_spikes ) );
    
    rec_spikes = [ net.rec_spikes, inf(2,num_spikes) ];
    num_rs = size(net.rec_spikes,2);
    
    Z = zeros(2, num_spikes, 'single');
    P = zeros(net.num_neurons+1, num_spikes, 'single');
    net.At = zeros(2,num_spikes, 'single');
    
    last_input_spikes = repmat(t, net.num_inputs, 1);
    last_output_spikes = repmat(t, net.num_neurons, 1);
    
    hX = net.hX;
    hZ = net.hZ;    
    
    W_exp = exp( double( net.W ) );
    V_n_exp = exp( -double( net.V ) );
    
    net.d_W =  zeros(net.num_neurons,net.num_inputs);
    net.d_V =  zeros(net.num_neurons,net.num_neurons);
    net.d_V0 = zeros(net.num_neurons,1);
    
    i = 1;
    l = 1;
    
    hX_all = zeros(net.num_inputs,num_spikes);
    hZ_all = zeros(net.num_neurons,num_spikes);
    
    if net.use_variance_tracking
        
        SW_new = net.SW;
        QW_new = net.QW;

        SV_new = net.SV;
        QV_new = net.QV;

        S0_new = net.S0;
        Q0_new = net.Q0;
        
        net.eta_W = net.eta*(QW_new-SW_new.^2)./(exp(-SW_new)+1);
        net.eta_V = net.eta*(QV_new-SV_new.^2)./(exp(-SV_new)+1);
        net.eta_0 = net.eta*(Q0_new-S0_new.^2)./(exp(-S0_new)+1);
    end
        
    try
        
        for j = 1:num_spikes

            t = spike_times(j);

            while (i < size(data.Xt,2)) && (t > data.Xt(2,i))

                n_id = data.Xt(1,i);
                sp_t = data.Xt(2,i);

                hX(n_id,1) = hX(n_id,1)*exp(-double(sp_t-last_input_spikes(n_id))/net.tau_x_r) + 1;
                hX(n_id,2) = hX(n_id,2)*exp(-double(sp_t-last_input_spikes(n_id))/net.tau_x_f) + 1;

                last_input_spikes(n_id) = sp_t;
                i = i+1;
            end
            
            while (l < size(rec_spikes,2)) && (t > rec_spikes(2,l))

                n_id = rec_spikes(1,l);
                sp_t = rec_spikes(2,l);
                
                hZ(n_id,1) = hZ(n_id,1)*exp(-double(sp_t-last_output_spikes(n_id))/net.tau_z_r) + 1;
                hZ(n_id,2) = hZ(n_id,2)*exp(-double(sp_t-last_output_spikes(n_id))/net.tau_z_f) + 1;

                last_output_spikes(n_id) = sp_t;
                l = l+1;
            end


            hZ(:,1) = hZ(:,1).*exp(-double(t-last_output_spikes)/net.tau_z_r);
            hZ(:,2) = hZ(:,2).*exp(-double(t-last_output_spikes)/net.tau_z_f);
            hX(:,1) = hX(:,1).*exp(-double(t-last_input_spikes)/net.tau_x_r);
            hX(:,2) = hX(:,2).*exp(-double(t-last_input_spikes)/net.tau_x_f);

            d_hX = diff(hX,1,2);
            d_hZ = diff(hZ,1,2);
            
            hX_all(:,j) = d_hX;
            hZ_all(:,j) = d_hZ;
            
            last_input_spikes(:) = t;
            last_output_spikes(:) = t;
            
            u_rf = net.w_rf*exp(-double(t-net.last_spike_t)/net.tau_rf);

            U = net.W*d_hX + net.V*d_hZ + u_rf;
            
            [P_t,A] = wta_softmax( U );

            [Z_i,k] = wta_draw_k(P_t);

            Z(:,j) = [k;t];
            P(:,j) = [P_t;t];
            
            rec_spikes(:,j+num_rs) = [k,t+net.rec_delay(k)];
            
            net.num_o(k) = net.num_o(k)+1;
            
            net.last_spike_t(k) = t;
            
            net.d_W(k,:) = net.d_W(k,:) + net.eta_W(k,:).*(d_hX'-W_exp(k,:))./max(net.eta_W(k,:),W_exp(k,:));
            
            net.d_V(:,k) = net.d_V(:,k) - net.eta_V(:,k);
            net.d_V(k,:) = net.d_V(k,:) + net.eta_V(k,:).*min(V_n_exp(k,:).*d_hZ',net.eta_V(k,:));

            net.At(:,j) = [A,t];
            
            if net.use_variance_tracking
                
                SW_new(k,:) = SW_new(k,:) + net.eta_W(k,:).*(net.W(k,:)+net.d_W(k,:)-SW_new(k,:));
                QW_new(k,:) = QW_new(k,:) + net.eta_W(k,:).*((net.W(k,:)+net.d_W(k,:)).^2-QW_new(k,:));
                
                SV_new(k,:) = SV_new(k,:) + net.eta_V(k,:).*((net.V(k,:)+net.d_V(k,:))-SV_new(k,:));
                QV_new(k,:) = QV_new(k,:) + net.eta_V(k,:).*((net.V(k,:)+net.d_V(k,:)).^2-QV_new(k,:));
            end
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
    
    if net.use_variance_tracking
        net.SW_new = SW_new;
        net.QW_new = QW_new;

        net.SV_new = SV_new;
        net.QV_new = QV_new;

        net.S0_new = S0_new;
        net.Q0_new = Q0_new;
    end
    
    net.rec_spikes = rec_spikes(:,l:end);
    net.rec_spikes(2,:) = net.rec_spikes(2,:)-t;
    net.last_spike_t = net.last_spike_t - t;
    
    net.R = sum( net.At(1,:) );
    
    net.hX = hX;
    net.hZ = hZ;
    net.hX_all = hX_all;
    net.hZ_all = hZ_all;

end
