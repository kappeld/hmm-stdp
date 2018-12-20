function [net,Z,P] = snn_sample_dthmm( net, data, ct )
% Draws a spike from sem network - discrete time spiking hmm sems realistic refractory mechanim.
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
%   tau_i_r         0.002        time constant of epsp window rise
%   tau_x_f         0.02         time constant of epsp window fall
%   tau_z_f         0.02         time constant of epsp window fall
%   tau_i_f         0.02         time constant of epsp window fall
%   tau_rf          0.005        time constant of refractory window
%   w_rf            -10          strength of refractroy window
%   w_i             4            strength of inhibitory spike
%   lambda          2000         network spike rate
%   rec_delay       0.005        delay on recurrent synapses
%   max_inhi_rate   1000         maximum rate of inhibitory neuron
%   sample_rate     2000         time discretisation, samples per second 
%   reset_psps      true         reset psps at beginning of sequence
%   use_inhibition  true         use correct inhibition for inference
%   use_variance_tracking  false
%
% @fields:
%   hX_all zeros   [num_inputs,sample_rate]    all feed forward psps
%   hZ_all zeros   [num_neurons,sample_rate]   all lateral psps
%   rf_all zeros   [num_neurons,sample_rate]   all feed forward psps
%   i_all  zeros   [1,sample_rate]             all inhibitory currents
%   At    szeros  [3,0]             log activity of network at spike times
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

    dt = 1/net.sample_rate;
    update_times = data.time(ct(1))+dt:dt:data.time(ct(end));
    
    num_spikes = 0;
    
    Z = zeros( 2, length(update_times), 'single' );
    P = zeros( net.num_neurons+1, length(update_times), 'single' );
    net.At = zeros( 2, length(update_times), 'single' );
    
    W_exp = exp( double( net.W ) );
    %W_n_exp = exp( -double( net.W ) );
    V_n_exp = exp( -double( net.V ) );
    V0_exp = exp( double( net.V0 ) );
    
    net.d_W =  zeros(net.num_neurons,net.num_inputs);
    net.d_V =  zeros(net.num_neurons,net.num_neurons);
    net.d_V0 = zeros(net.num_neurons,1);
    
    X_psp = exp(-(0:dt:(10*net.tau_x_f))/net.tau_x_f) - exp(-(0:dt:(10*net.tau_x_f))/net.tau_x_r);
    Z_psp = exp(-(0:dt:(10*net.tau_z_f))/net.tau_z_f) - exp(-(0:dt:(10*net.tau_z_f))/net.tau_z_r);
    I_psp = exp(-(0:dt:(10*net.tau_i_f))/net.tau_i_f) - exp(-(0:dt:(10*net.tau_i_f))/net.tau_i_r);
    R_psp = exp(-(0:dt:(10*net.tau_rf))/net.tau_rf);

    delay_samples = round( net.rec_delay*net.sample_rate );
    
    num_prev_samples = ceil( net.sample_rate*( 10*max( [net.tau_x_f,net.tau_z_f,net.tau_rf] ) ) ) + ...
                       delay_samples + 1;
                                           
    hX_all = zeros( net.num_inputs, length(update_times) + num_prev_samples );
    hZ_all = zeros( net.num_neurons, length(update_times) + num_prev_samples );
    rf_all = zeros( net.num_neurons, length(update_times) + num_prev_samples );
    i_all = zeros( 1, length(update_times) + num_prev_samples );
    
    if ~net.reset_psps
        hX_all(:,1:num_prev_samples) = net.hX_all(:,(end-num_prev_samples+1):end);
        hZ_all(:,1:num_prev_samples) = net.hZ_all(:,(end-num_prev_samples+1):end);
        rf_all(:,1:num_prev_samples) = net.rf_all(:,(end-num_prev_samples+1):end);
        i_all(:,1:num_prev_samples) = net.i_all(:,(end-num_prev_samples+1):end);
    end
    
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
    
    i = 1;
    
    try        
        for j = 1:length(update_times)

            t = update_times(j);

            while (i < size(data.Xt,2)) && (t > data.Xt(2,i))
                n_id = data.Xt(1,i);
                i = i+1;
                hX_all(n_id,j+(1:length(X_psp))) = hX_all(n_id,j+(1:length(X_psp))) + X_psp;
                
                %net.d_W(:,n_id) = net.d_W(:,n_id) - net.eta_W(:,n_id);
            end
            
            d_hX = hX_all(:,j);
            d_hZ = hZ_all(:,j);

            %U = (1/60)*net.W*d_hX + (1/4)*net.V*d_hZ + net.w_rf*rf_all(:,j) - net.V0; % + net.w_i*i_all(j);
            U = net.W*d_hX + net.V*d_hZ + net.w_rf*rf_all(:,j);
            
            [P_t,A] = wta_softmax( U );
            
            if ~net.use_inhibition
                P_t = exp( U );
            end
            
            P(:,j) = [P_t;t];
            
            Z_i = exprnd(1./(net.lambda*P_t)) < dt;
            
            %Z_i = exprnd(1./(net.lambda*(min(1,exp(U))))) < dt;
            
            %U_i = sum( d_hZ );
            
            net.At(:,j) = [A,t];
            
            %if exprnd(1./(net.max_inhi_rate*(min(1,exp(U_i))))) < dt
            %    i_all(:,j+(1:length(X_psp))) = i_all(:,j+(1:length(I_psp))) + I_psp;
            %end

            for k = find(Z_i)'
            
                net.num_o(k) = net.num_o(k)+1;

                num_spikes = num_spikes+1;
                Z(:,num_spikes) = [k,t];
                
                hZ_all(k,delay_samples+j+(1:length(Z_psp))) = ...
                    hZ_all(k,delay_samples+j+(1:length(Z_psp))) + Z_psp;
                
                rf_all(k,j+(1:length(R_psp))) = rf_all(k,j+(1:length(R_psp))) + R_psp;

                net.d_W(k,:) = net.d_W(k,:) + net.eta_W(k,:).*(d_hX'-W_exp(k,:))./max(net.eta_W(k,:),W_exp(k,:));
                %net.d_W(k,:) = net.d_W(k,:) + net.eta_W(k,:).*min(W_n_exp(k,:).*d_hX',net.eta_W(k,:));
                
                net.d_V(:,k) = net.d_V(:,k) - net.eta_V(:,k);
                net.d_V(k,:) = net.d_V(k,:) + net.eta_V(k,:).*min(V_n_exp(k,:).*d_hZ',net.eta_V(k,:));
                
                if net.use_variance_tracking

                    SW_new(k,:) = SW_new(k,:) + net.eta_W(k,:).*(net.W(k,:)+net.d_W(k,:)-SW_new(k,:));
                    QW_new(k,:) = QW_new(k,:) + net.eta_W(k,:).*((net.W(k,:)+net.d_W(k,:)).^2-QW_new(k,:));

                    SV_new(k,:) = SV_new(k,:) + net.eta_V(k,:).*((net.V(k,:)+net.d_V(k,:))-SV_new(k,:));
                    QV_new(k,:) = QV_new(k,:) + net.eta_V(k,:).*((net.V(k,:)+net.d_V(k,:)).^2-QV_new(k,:));
                end
            end
            
            if any( Z_i )
                
                net.d_V0 = net.d_V0 + net.eta_0.*(Z_i-V0_exp)./max(net.eta_0,V0_exp);

                if net.use_variance_tracking

                    S0_new = S0_new + net.eta_0.*((net.V0+net.d_V0)-S0_new);
                    Q0_new = Q0_new + net.eta_0.*((net.V0+net.d_V0).^2-Q0_new);
                end
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
    
    Z = Z(:,1:num_spikes);
    
    net.R = sum( net.At(1,:) );
    
    net.hX_all = hX_all;
    net.hZ_all = hZ_all;
    net.rf_all = rf_all;
    net.i_all = i_all;

end
