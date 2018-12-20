function [net,Z,P] = snn_sample_ct_normal( net, data, ct )
% Draws a spike from sem network - continuous time spiking hmm sems.
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
%   tau_x         0.02         time constant of epsp window
%   tau_z         0.02         time constant of epsp window
%   lambda        2000         network spike rate
%   refrac_time   0.002        absolute refractory time
%
% @fields:
%   At    szeros  [2,0]             log activity of network at spike times
%   hX    zeros   [num_inputs,1]    feedforward PSPs
%   hZ    zeros   [num_neurons,1]   recurrent PSPs
%   R     zeros   [1,1]             current importance weight
%   

    t = data.time(ct(1));
    i = 1;
    
    %if isfield( net, 'Delta' )
    %    Delta = net.Delta;
    %    num_spikes = length(Delta);
    %else
        time_range = data.time(ct(end))-data.time(ct(1));    
        num_spikes = poissrnd( double( net.lambda*time_range ) );
        Delta = diff( [0, sort( time_range*rand( 1, num_spikes ) ) ] );
        net.Delta = Delta;
    %end
    
    Z = zeros(2, num_spikes, 'single');
    P = zeros(net.num_neurons+1, num_spikes, 'single');
    net.At = zeros(2,num_spikes, 'single');
    
    last_input_spikes = repmat(t, net.num_inputs, 1);
    
    hX = net.hX;
    hZ = net.hZ;    
    eta = net.eta;
    
    W = double( net.W );
    V = double( net.V );
    V0_exp = exp( double( net.V0 ) );
    
    net.d_W =  zeros(net.num_neurons,net.num_inputs);
    net.d_V =  zeros(net.num_neurons,net.num_neurons);
    net.d_V0 = zeros(net.num_neurons,1);
    
    %if ~isfield( net, 'hX_all' );
    %    net.hX_all = [];
    %    net.hZ_all = [];
    %end
    
    try
        
        for j = 1:length(Delta)

            t = t+Delta(j);

            while (i < size(data.Xt,2)) && (t > data.Xt(2,i))

                n_id = data.Xt(1,i);
                sp_t = data.Xt(2,i);

                hX(n_id) = hX(n_id)*exp(-(sp_t-last_input_spikes(n_id))/net.tau_x) + 1;

                last_input_spikes(n_id) = sp_t;

                i = i+1;
            end

            hZ = hZ*exp(-Delta(j)/net.tau_z);
            hX = hX.*exp(-(t-last_input_spikes)/net.tau_x);

            last_input_spikes(:) = t;

            U = net.W*hX + net.V*hZ + net.V0;
            [P_t,A] = wta_softmax( U );

            [Z_i,k] = wta_draw_k(P_t);

            Z(:,j) = [k;t];
            P(:,j) = [P_t;t];

            net.d_W(k,:) = net.d_W(k,:) + hX' - W(k,:);
            net.d_V(k,:) = net.d_V(k,:) + hZ' - V(k,:);
            net.d_V0     = net.d_V0     + (Z_i-V0_exp)./max(eta,V0_exp);

            net.At(:,j) = [A,t];
            
            %net.hX_all = [net.hX_all,hX];
            %net.hZ_all = [net.hZ_all,hZ];

            hZ(k) = hZ(k) + 1;
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
    
    net.R = sum( net.At(1,:) );    
    net.hX = hX;
    net.hZ = hZ;

end
