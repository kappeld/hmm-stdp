function [net,Z,P] = snn_sample_dtrf( net, data, ct )
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
%   tau_x_r         0.002        time constant of epsp window rise [s]
%   tau_z_r         0.002        time constant of epsp window rise [s]
%   tau_x_f         0.02         time constant of epsp window fall [s]
%   tau_z_f         0.02         time constant of epsp window fall [s]
%   tau_rf          0.005        time constant of refractory window [s]
%   tau_I_r         0.002        time constant of Inhibitory window rise [s]
%   tau_I_f         0.010        time constant of Inhibitory window fall [s]
%   w_rf            -10          strength of refractroy window
%   lambda          2000         network spike rate [1/s]
%   rec_delay       0.005        delay on recurrent synapses [s]
%   sample_rate     2000         time discretisation [1/s]
%   reset_psps      true         reset psps at beginning of sequence
%   I_int_time      0.01         integration time of inhibitory current [s]
%   I_delay         0.005        delay of inhibitory current [s]
%   I_refrac_time   0.005        refractory time of inhibitory neuron [s]
%   hom_ts          0.001        track speed of locale homeostatsis
%   eta             0.005        learn rate
%   learn_delay     0.0          delay of weight consolidation [s]
%
% @fields:
%   hX_all zeros      [num_inputs,sample_rate]    all feed forward psps
%   hZ_all zeros      [num_neurons,sample_rate]   all lateral psps
%   rf_all zeros      [num_neurons,sample_rate]   all feed forward psps
%   I_all zeros       [1,sample_rate]             all inhibitory psps
%   R_all zeros       [1,sample_rate]             all importance weights
%   At    szeros      [3,0]                       log activity of network at spike times
%   R     zeros       [1,1]                       current importance weight
%   I_s   zeros       [num_neurons,1]             current slow inhibition
%   I_f   zeros       [1,1]                       current fast inhibition
%   U_w_offs   zeros  [num_neurons,1]             offset of recurrent connections
%   U_v_offs   zeros  [num_neurons,1]             offset of feedforward connections
%   U_x_offs   ones   [1,1]                       offset of recurrent connections
%   U_z_offs   ones   [1,1]                       offset of feedforward connections
%   I_ref zeros       [1,1]                       refractory counter of inhibition
%   num_o zeros       [num_neurons,1]             number of output spikes per neuron
%

    dt = 1/net.sample_rate;
    update_times = data.time(ct(1))+dt:dt:data.time(ct(end));
    
    num_spikes = 0;
    
    Z = zeros( 2, length(update_times), 'single' );
    P = zeros( net.num_neurons+1, length(update_times), 'single' );
    net.At = zeros( 5, length(update_times), 'single' );
    
    W_exp = exp( double( net.W ) );
    V_exp = exp( double( net.V ) );
    V0_exp = exp( double( net.V0 ) );
    
    num_learn_samples = ceil( net.sample_rate*net.learn_delay );
    
    if ( size(net.d_W,3) > num_learn_samples )
        W_tmp = net.d_W(:,:,(end-num_learn_samples+1):end);
        net.d_W = zeros(net.num_neurons,net.num_inputs,length(update_times)+num_learn_samples);
        net.d_W(:,:,1:num_learn_samples) = W_tmp;
        
        V_tmp = net.d_V(:,:,(end-num_learn_samples+1):end);
        net.d_V = zeros(net.num_neurons,net.num_neurons,length(update_times)+num_learn_samples);
        net.d_V(:,:,1:num_learn_samples) = V_tmp;
        
        net.d_V0 = [ net.d_V0(:,(end-num_learn_samples+1):end), zeros(net.num_neurons,length(update_times)) ];
    else
        net.d_W = zeros(net.num_neurons,net.num_inputs,length(update_times)+num_learn_samples);
        net.d_V = zeros(net.num_neurons,net.num_neurons,length(update_times)+num_learn_samples);
        net.d_V0 = zeros(net.num_neurons,1,length(update_times)+num_learn_samples);
    end
    
    X_psp = exp(-(0:dt:(10*net.tau_x_f))/net.tau_x_f) - exp(-(0:dt:(10*net.tau_x_f))/net.tau_x_r);
    Z_psp = exp(-(0:dt:(10*net.tau_z_f))/net.tau_z_f) - exp(-(0:dt:(10*net.tau_z_f))/net.tau_z_r);
    R_psp = exp(-(0:dt:(10*net.tau_rf))/net.tau_rf);
    I_psp = 5*( exp(-(0:dt:(10*net.tau_I_f))/net.tau_I_f) - exp(-(0:dt:(10*net.tau_I_f))/net.tau_I_r) );
    
    num_prev_samples = ceil( net.sample_rate*( 10*max( [net.tau_x_f, net.tau_z_f, net.tau_rf, net.tau_I_f] ) + ...
                                               max( net.rec_delay ) ) ) + 5;
                                           
    I_delay_samples = ceil( net.I_delay*net.sample_rate );
                                           
    hX_all = zeros( net.num_inputs, length(update_times) + num_prev_samples );
    hZ_all = zeros( net.num_neurons, length(update_times) + num_prev_samples );
    rf_all = zeros( net.num_neurons, length(update_times) + num_prev_samples );
    I_all = zeros( 1, length(update_times) + num_prev_samples );
    R_all = zeros( 1, length(update_times) );
    
    if ~net.reset_psps
        hX_all(:,1:num_prev_samples) = net.hX_all(:,(end-num_prev_samples+1):end);
        hZ_all(:,1:num_prev_samples) = net.hZ_all(:,(end-num_prev_samples+1):end);
        rf_all(:,1:num_prev_samples) = net.rf_all(:,(end-num_prev_samples+1):end);
        I_all(:,1:num_prev_samples) = net.I_all(:,(end-num_prev_samples+1):end);
    end
    
    delay_samples = round( net.rec_delay*net.sample_rate );
    I_ref_samples = max( 0, round( net.I_refrac_time*net.sample_rate - 1 ) );

    if ~isnan( net.self_inhibition )
        net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
    end
    
    hom_ts = net.hom_ts;
    
    i = 1;
    
    try        
        for j = 1:length(update_times)

            t = update_times(j);

            while (i < size(data.Xt,2)) && (t > data.Xt(2,i))
                n_id = data.Xt(1,i);
                i = i+1;
                hX_all(n_id,j+(1:length(X_psp))) = hX_all(n_id,j+(1:length(X_psp))) + X_psp;
            end
            
            d_hX = hX_all(:,j);
            d_hZ = hZ_all(:,j);
            u_rf = rf_all(:,j);
            
            % keep mean of W and V part of membrane potential at zero
            W_offs = repmat( net.U_w_offs/net.U_x_offs, 1, net.num_inputs );
            V_offs = repmat( net.U_v_offs/net.U_z_offs + 1, 1, net.num_neurons );            
            U_w = (net.W - W_offs)*d_hX;
            U_v = (net.V - V_offs)*d_hZ + net.w_rf*u_rf + net.V0;
            
            U = U_w + U_v - I_all(j) - net.I_s;
            
            % assure that contirbution of w and v for each  neuron is 0 on
            % average. This should give a good homeostatic effect.
            %net.U_w_offs = (1-hom_ts)*net.U_w_offs + hom_ts*net.W*d_hX;
            %net.U_v_offs = (1-hom_ts)*net.U_v_offs + hom_ts*net.V*d_hZ;
            %net.U_x_offs = (1-hom_ts)*net.U_x_offs + hom_ts*sum(d_hX);
            %net.U_z_offs = (1-hom_ts)*net.U_z_offs + hom_ts*sum(d_hZ);

            P_t = exp(U);
            
            A_t = log( sum( P_t ) );

            P(:,j) = [P_t;t];
            
            net.At(:,j) = [A_t;I_all(j);mean(net.I_s);exp(net.R-net.r_mean);t];

            Z_i = exprnd(1./(net.lambda*P_t)) < dt;
            
            net.R = net.R*0.99 + 0.01*max(-1000,log(I_all(j)))/50;
            
            R_all(j) = net.R;
            
            % globale homeostasis
            net.I_s = net.I_s + 0.005*(sum(Z_i) - 200/net.sample_rate);
            
            % locale homeostasis
            %net.I_s = net.I_s + 0.01*(Z_i - 50/net.sample_rate);
            
            if (net.I_ref > 0)
                net.I_ref = net.I_ref-1;
            elseif any( Z_i )
                I_all(I_delay_samples+j+(1:length(I_psp))) = ...
                    I_all(I_delay_samples+j+(1:length(I_psp))) + I_psp;
                    
                net.I_ref = I_ref_samples;
            end
            
            w_j = num_learn_samples + j;
            
            for k = find(Z_i)'
                
                net.num_o(k) = net.num_o(k)+1;

                num_spikes = num_spikes+1;
                Z(:,num_spikes) = [k,t];
                
                hZ_all(k,delay_samples+j+(1:length(Z_psp))) = ...
                    hZ_all(k,delay_samples+j+(1:length(Z_psp))) + Z_psp;
                
                rf_all(k,j+(1:length(R_psp))) = rf_all(k,j+(1:length(R_psp))) + R_psp;

                net.d_W(k,:,w_j) = net.d_W(k,:,w_j) + net.eta.*(d_hX'-W_exp(k,:))./max(net.eta,W_exp(k,:));
                net.d_V(k,:,w_j) = net.d_V(k,:,w_j) + net.eta.*(d_hZ'-V_exp(k,:))./max(net.eta,V_exp(k,:));
            end
            
            if any( Z_i )
                net.d_V0(:,w_j) = net.d_V0(:,w_j) + net.eta.*(Z_i-V0_exp)./max(net.eta,V0_exp);
            end
            
            % delayed update of network weights
            net.r_mean = max( net.r_mean, R_all(j) );

            p_accept = 1; %exp( double(R_all(j)) - net.r_mean );

            net.W = net.W + p_accept*net.d_W(:,:,j);
            net.V = net.V + p_accept*net.d_V(:,:,j);
            net.V0 = net.V0 + p_accept*net.d_V0(:,j);

            if ~isnan( net.self_inhibition )
                net.V(eye(net.num_neurons) > 0) = net.self_inhibition;
            end            
        end
    catch
        fprintf('There has been an error, while sampling!\nExcluding run from training\n');
        file_name = sprintf('/tmp/error_data_%u_%04u.mat', net.iteration, round(9999*rand()) );
        the_error = lasterror();
        fprintf( '\n%s\n', the_error.message );
        fprintf( '  in %s on line %i\n\n', the_error.stack(1).name, the_error.stack(1).line );
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
    
    net.hX_all = hX_all;
    net.hZ_all = hZ_all;
    net.rf_all = rf_all;
    net.I_all = I_all;
    net.R_all = R_all;
end
