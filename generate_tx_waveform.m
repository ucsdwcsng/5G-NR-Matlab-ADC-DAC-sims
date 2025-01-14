clear all
%% Important params
num_frames = 1;
SCS = 30;
grid_size = 273; % num RBs, 273 corr

num_layers = 1;
num_tx_ant = 1;
num_rx_ant = 1;

num_slots_to_gen = 10;
num_zeros_to_append = 10000;
num_zeros_to_prepend = 10000; % zeros

[tx_samples_ideal,sim_params,waveform_params] = get_tx_iq_samples(num_slots_to_gen,SCS,grid_size,num_layers,num_tx_ant,num_rx_ant);
