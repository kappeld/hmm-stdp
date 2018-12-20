function net = snn_train_rs_lin( net, data, ct )
% Spiking EM learning - continuous time spiking hmm sems - rejection sampling.
%
% net = snn_train_ct( net, data, ct )
%
%
% @parameters:
%   sample_method       'ct'    default sample method
%   performance_method  'ct'    default performance method
%   self_inhibition     0       fix value of self-recurrent weights
%   use_iw              true    use importance weight
%   resample            true    resample path if rejected
%   p_required_fields   [R_all,d_W,d_V,d_V0]
%   p_requirements      [batch_mode]
%
% @fields:
%   W        rand[-0.01]    [num_neurons,num_inputs]  feedforward network weights
%   V        rand[-0.01]    [num_neurons,num_neurons] recurrent network weights
%   V0       rand[-0.01]    [num_neurons,1]           network bias
%   d_W      zeros          [num_neurons,num_inputs]  feedforward weight change
%   d_V      zeros          [num_neurons,num_neurons] recurrent weight change
%   d_V0     zeros          [num_neurons,1]           bias change
%   r_mean   const[-inf]    [1,1]                     current estimate of iw normalisation
%
%
% David Kappel
% 24.05.2011
%
%

    net.iteration = net.iteration + 1;
    net.got_sample = true;
end
