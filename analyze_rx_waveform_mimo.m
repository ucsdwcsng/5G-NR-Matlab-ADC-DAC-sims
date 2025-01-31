clear all
close all
%% Important params
rng(2)
num_frames = 1;
SCS = 30;
grid_size = 273; % num RBs, 273 corr 100 mhz
% grid_size = 72; % num RBs, 273 corr
% grid_size = 52; % num RBs, 273 corr 30 mhz

num_layers = 2;
num_tx_ant = 2;
num_rx_ant = 2;

num_slots_to_gen = 1;
% num_zeros_to_append = 10000;
% num_zeros_to_prepend = 10000; % zeros

get_tx_params = true;
if(get_tx_params)
    [tx_samples_ideal,sim_params,waveform_params] = get_tx_iq_samples_mimo(num_slots_to_gen,SCS,grid_size,num_layers,num_tx_ant,num_rx_ant);
    num_repmats = 1;
    tx_samples_repeated_with_zeros = repmat(tx_samples_ideal,[1,num_repmats]);
    tx_samples_repeated_with_zeros = padarray(tx_samples_repeated_with_zeros.',[100,0],0,'post'); % little bit extra
end


% pspectrum(tx_samples_ideal,waveform_params.SampleRate)
% waveform_params.SampleRate
%  write_complex_binary(tx_samples_ideal,"5g_ofdm_tx_norm.iq")
%% Load RX waveform from aux hardwre (Code goes here, or manual load)

nominal_SNR = 25;
SNR = 10^((nominal_SNR+16)/10);
N0 = (1/sqrt(SNR));
% nPowerPerRE = N0^2*double(waveform_params.Nfft);
noise = N0*randn(size(tx_samples_repeated_with_zeros),"like",tx_samples_repeated_with_zeros);
% noise_power = 20*log10(mean(abs(noise)))
% sig_power = 20*log10(mean(abs(tx_samples_repeated_with_zeros)))

phase_offset_between_antennas = 60; % degrees

tx_samples_with_artf_channel = tx_samples_repeated_with_zeros;
% tx_samples_with_artf_channel(:,1) = tx_samples_repeated_with_zeros(:,1)*exp(1j*(pi/180)*phase_offset_between_antennas);
% tx_samples_with_artf_channel(:,2) = tx_samples_repeated_with_zeros(:,2);
% tx_samples_with_artf_channel(:,3) = tx_samples_repeated_with_zeros(:,3)*exp(1j*(pi/180)*phase_offset_between_antennas);
% tx_samples_with_artf_channel(:,4) = tx_samples_repeated_with_zeros(:,4);

rxWaveform = tx_samples_with_artf_channel+noise;
enable_artificial_cfo = false;
if(enable_artificial_cfo)
    f_offset_khz = 1e3;
    % delta_f = 0.001*f_offset_khz;
    % delta_f=0;
    Fs = waveform_params.SampleRate;
    % f_offset = (f_offset_khz+delta_f*(rand(1,length(rxWaveform))-0.5))/Fs; % can go as high as 1e-4
    f_offset = (f_offset_khz)/Fs; % can go as high as 1e-4
    rxWaveform = rxWaveform.*(exp(-1i*2*pi*(f_offset).*(0:length(rxWaveform)-1)).');
end
IQImbalanceON = false;
if IQImbalanceON
    amplitudeImbalance = 0.2;
    phaseImbalance = 0.5;
    rxWaveform = iqimbal(rxWaveform,amplitudeImbalance,phaseImbalance);
end

rxWaveform_fin = rxWaveform;


[mean_SNR_per_layer,snr_diff,sim_params] = analyze_rx_iq_samps(rxWaveform_fin,sim_params,waveform_params); 

mean_SNR_per_layer
snr_diff