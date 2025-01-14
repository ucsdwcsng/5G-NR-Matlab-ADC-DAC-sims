clear all
%% Important params
rng(1)
num_frames = 1;
SCS = 30;
grid_size = 273; % num RBs, 273 corr
% grid_size = 52; % num RBs, 273 corr

num_layers = 1;
num_tx_ant = 1;
num_rx_ant = 1;

write_tx_to_file_rfsoc = 1; %enable this parameter to write the data into mem file for RF-SoC transmission.

num_slots_to_gen = 1;
% num_zeros_to_append = 10000;
% num_zeros_to_prepend = 10000; % zeros

get_tx_params = true;
if(get_tx_params)
    [tx_samples_ideal,sim_params,waveform_params] = get_tx_iq_samples(num_slots_to_gen,SCS,grid_size,num_layers,num_tx_ant,num_rx_ant);
    if(~write_tx_to_file_rfsoc)
        num_repmats = 10;
        tx_samples_ideal = repmat(tx_samples_ideal,[1,num_repmats]);
    end
end
% pspectrum(tx_samples_ideal,waveform_params.SampleRate)
% waveform_params.SampleRate

if(write_tx_to_file_rfsoc == 1)
    sf = 2^14;
    tx_samples = resample(tx_samples_ideal.',4,1); % Data is upsampled by 4X for RF-SoC Transmission.
    tx_samples_hex = cell(size(tx_samples,1), 1);
    fileID = fopen('tx_samples_hex_5G.mem', 'w');
    
    for k = 1:size(tx_samples,1)
        iq = tx_samples(k);
        tx_samples_hex{k} = [dec2hex(round(imag(iq*sf)),4),dec2hex(round(real(iq*sf)),4)];
    end
    
    for k=1:4:size(tx_samples,1)
        fprintf(fileID,'%s%s%s%s\n',tx_samples_hex{k+3},tx_samples_hex{k+2},tx_samples_hex{k+1},tx_samples_hex{k});
    end

    fclose(fileID);

end

%  write_complex_binary(tx_samples_ideal,"5g_ofdm_tx_norm.iq")
%% Load RX waveform from aux hardwre (Code goes here, or manual load)
sim = true;
% store_to_file = true;
if(sim)
    tx_samples_ideal = [tx_samples_ideal zeros(1,100)]; % leeway
    nominal_SNR = 34;
    SNR = 10^((nominal_SNR-17)/10);
    N0 = (1/sqrt(num_rx_ant*double(waveform_params.Nfft)*(SNR)));
    nPowerPerRE = N0^2*double(waveform_params.Nfft);
    noise = N0*randn(size(tx_samples_ideal),"like",tx_samples_ideal);
    rxWaveform = tx_samples_ideal+noise;
    enable_artificial_cfo = true;
    if(enable_artificial_cfo)
        f_offset_khz = 1e3;
        % delta_f = 0.001*f_offset_khz;
        % delta_f=0;
        Fs = waveform_params.SampleRate;
        % f_offset = (f_offset_khz+delta_f*(rand(1,length(rxWaveform))-0.5))/Fs; % can go as high as 1e-4
        f_offset = (f_offset_khz)/Fs; % can go as high as 1e-4
        rxWaveform = rxWaveform.*exp(-1i*2*pi*(f_offset).*(0:length(rxWaveform)-1));
    end
    IQImbalanceON = false;
    if IQImbalanceON
        amplitudeImbalance = 0.2;
        phaseImbalance = 0.5;
        rxWaveform = iqimbal(rxWaveform,amplitudeImbalance,phaseImbalance);
    end

    % save("./iq-files/matlab_data/rxWaveform.mat","rxWaveform")
else
    % captured_data=load("./Data/01_12_24/capture_rx0.txt");
    replay_hwr_trace = true;

    if(replay_hwr_trace)
        % captured_data=load("./Data/10_12_24/5g_ofdm_norm_40gain_4M_12x.txt");
        captured_data = load("./Data/11_12_24/40gain_4M-8x.txt");
        % captured_data=load("capture_rx0.txt");
        disp("Captured data")
        captured_data=captured_data(1:1:end,1:2);
        % maxi=max(max(captured_data));
        % mini=min(min(captured_data));
        % N=log2(maxi-mini);
        N=16;
    
        i_samp= bit2vol(captured_data(:,1),N);
        q_samp= bit2vol(captured_data(:,2),N);
    
        rx_samples=i_samp+1i*q_samp;
        disp("Complex data reconstructed")
        rxWaveform = resample(rx_samples,1,8);
        plot(abs(rxWaveform))
    else
        load("./iq-files/matlab_data/rxWaveform.mat")
    end
end

%% Analyze the rx waveform
rxWaveform = resample(conj(iq_adc*exp(1j*3*(pi/2))),1,4);
% rxWaveform_trun = rxWaveform(1:1:40000);
[mean_SNR_per_slot,sim_params] = analyze_rx_iq_samps(rxWaveform.',sim_params,waveform_params); 