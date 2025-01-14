function [output_snr,output_snr_diff,diag_snr,...
    diag_snr_diff,ref_snr,ref_snr_diff] ...
    = output_snr_recovered(input_snr,seed,interp_fac)
    %% Sim params
    tic
    rng(seed)
    num_frames = 1;
    SCS = 30;
    grid_size = 273; % num RBs, 273 corr 100 mhz
    SNR = input_snr;
    % grid_size = 72; % num RBs, 273 corr
    % grid_size = 52; % num RBs, 273 corr 30 mhz
    
    num_layers = 4;
    num_tx_ant = 4;
    num_rx_ant = 4;
    
    num_slots_to_gen = 1;
    % num_zeros_to_append = 10000;
    % num_zeros_to_prepend = 10000; % zeros
    
    [tx_samples_ideal,sim_params,waveform_params] = get_tx_iq_samples_mimo(num_slots_to_gen,SCS,grid_size,num_layers,num_tx_ant,num_rx_ant);
    num_repmats = 1;
    tx_samples_repeated_with_zeros = repmat(tx_samples_ideal,[1,num_repmats]);
    txWaveform = padarray(tx_samples_repeated_with_zeros.',[100,0],0,'post'); % little bit extra

    
    %%
    analog_factor = 160;
    adc_oversamp_factor_trad = 4;
    adc_oversamp_factor = 16;
    
    rxWaveform_analog_afterLNA_nonoise = resample(txWaveform,analog_factor,1);
    noise_sigma = sqrt(0.5*10^(-(SNR+10)/10)); % +6 normalizes
    
    rx_waveform_size = size(rxWaveform_analog_afterLNA_nonoise);
    complex_noise =  noise_sigma*(randn(rx_waveform_size)+1j*randn(rx_waveform_size));
    rxWaveform_analog_afterLNA = rxWaveform_analog_afterLNA_nonoise+complex_noise;
    
    rxWaveform_ADC_Samples_raw = rxWaveform_analog_afterLNA(1:(analog_factor/adc_oversamp_factor):end,:); % sample and hold ADC
    rxWaveform_ADC_Samples_trad = rxWaveform_analog_afterLNA(1:(analog_factor/adc_oversamp_factor_trad):end,:); % sample and hold ADC
    
    rxWaveform_ADC_Samples_combined = zeros([size(rxWaveform_ADC_Samples_raw,1),1]);
    
    switch_vec = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4];
    
    for i=1:1:16
        rxWaveform_ADC_Samples_combined(i:16:end,1) = rxWaveform_ADC_Samples_raw(i:16:end,switch_vec(i));
    end
    rxWaveform_ADC_decim = lowpass(rxWaveform_ADC_Samples_combined,interp_fac);
    % rxWaveform_ADC_Samples_trad_decim = lowpass(rxWaveform_ADC_Samples_trad,0.25);
    rxWaveform_ADC1x = resample(rxWaveform_ADC_Samples_trad,1,adc_oversamp_factor_trad);
    % rxWaveform_ADC1x = rxWaveform_ADC_Samples_trad_decim(1:4:end,:);

    % rxWaveform_ADC_decim = rxWaveform_ADC_Samples_raw;
    % rxWaveform_DBF  = rxWaveform_ADC_decim(1:4:end,:);
    rxWaveform_GMO = zeros([4*size(txWaveform,1),4]);
    rxWaveform_GMO(1:4:end,1) = rxWaveform_ADC_decim(1:16:end,1);
    rxWaveform_GMO(2:4:end,1) = rxWaveform_ADC_decim(2:16:end,1);
    rxWaveform_GMO(3:4:end,1) = rxWaveform_ADC_decim(3:16:end,1);
    rxWaveform_GMO(4:4:end,1) = rxWaveform_ADC_decim(4:16:end,1);
    
    rxWaveform_GMO(1:4:end,2) = rxWaveform_ADC_decim(5:16:end,1);
    rxWaveform_GMO(2:4:end,2) = rxWaveform_ADC_decim(6:16:end,1);
    rxWaveform_GMO(3:4:end,2) = rxWaveform_ADC_decim(7:16:end,1);
    rxWaveform_GMO(4:4:end,2) = rxWaveform_ADC_decim(8:16:end,1);
    
    rxWaveform_GMO(1:4:end,3) = rxWaveform_ADC_decim(9:16:end,1);
    rxWaveform_GMO(2:4:end,3) = rxWaveform_ADC_decim(10:16:end,1);
    rxWaveform_GMO(3:4:end,3) = rxWaveform_ADC_decim(11:16:end,1);
    rxWaveform_GMO(4:4:end,3) = rxWaveform_ADC_decim(12:16:end,1);
    
    rxWaveform_GMO(1:4:end,4) = rxWaveform_ADC_decim(13:16:end,1);
    rxWaveform_GMO(2:4:end,4) = rxWaveform_ADC_decim(14:16:end,1);
    rxWaveform_GMO(3:4:end,4) = rxWaveform_ADC_decim(15:16:end,1);
    rxWaveform_GMO(4:4:end,4) = rxWaveform_ADC_decim(16:16:end,1);
    
    rxWaveform_ADC4x = resample(rxWaveform_GMO,1,4);
    % 
    % rxWaveform_ADC_Samples_after_switching = zeros([1,numel(rxWaveform_ADC_Samples_raw)]);
    % rxWaveform_ADC_Samples_after_switching
    % rxWaveform_analog_afterLNA(1:(analog_factor/adc_oversamp_factor):end,:); % sample and hold ADC
    
    % rxWaveform = resample(rxWaveform_ADC_Samples_raw,1,adc_oversamp_factor);
    % rxWaveform = resample(rxWaveform_analog_afterLNA,1,analog_factor);
    
    plot_filters = false;
    if(plot_filters)
        % analog signal
        Fs = carrier_bw;
        analog_freq_vec = linspace(-analog_factor*(Fs/2),analog_factor*(Fs/2),numel(rxWaveform_analog_afterLNA(:,1))); 
        adc_freq_vec = linspace(-adc_oversamp_factor*(Fs/2),adc_oversamp_factor*(Fs/2),numel(rxWaveform_ADC_Samples_raw(:,1))); 
        bb_freq_vec = linspace(-(Fs/2),(Fs/2),numel(rxWaveform_DBF(:,1))); 
    
        rxWaveform_analog_afterLNA_freq = fftshift(fft(rxWaveform_analog_afterLNA(:,1)));
        rxWaveform_ADC_Samples_raw_freq = fftshift(fft(rxWaveform_ADC_Samples_raw(:,1)))*(analog_factor/adc_oversamp_factor);
        rxWaveform_ADC_decim_freq = fftshift(fft(rxWaveform_ADC_decim(:,1)))*(analog_factor/adc_oversamp_factor);
        rxWaveform_freq = fftshift(fft(rxWaveform_DBF(:,1)))*analog_factor; % normalize
        figure;
        plot(bb_freq_vec,20*log10(abs(rxWaveform_freq)));
        hold on
        plot(adc_freq_vec,20*log10(abs(rxWaveform_ADC_decim_freq)));
        plot(adc_freq_vec,20*log10(abs(rxWaveform_ADC_Samples_raw_freq)));
        plot(analog_freq_vec,20*log10(abs(rxWaveform_analog_afterLNA_freq)));
        
        legend(fliplr(["After LNA","ADC Raw","After digital decim filter","Downsampled"]))
    end
    %%
    % rxWaveform = zeros(size(txWaveform));
    % rxWaveform(:,1) = txWaveform(:,1)+0.1*exp(1j*pi/2)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
    % rxWaveform(:,2) = txWaveform(:,2)+0.1*exp(1j*pi/18)*txWaveform(:,1)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4);
    % rxWaveform(:,3) = txWaveform(:,3)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,1)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
    % rxWaveform(:,4) = txWaveform(:,4)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,1); 
    % 
    % rxWaveform(:,1) = txWaveform(:,1); 
    % rxWaveform(:,2) = txWaveform(:,2);
    % rxWaveform(:,3) = txWaveform(:,3); 
    % rxWaveform(:,4) = txWaveform(:,4); 
    
    plot_spectrum = false;
    if(plot_spectrum)
        test_wf_time = rxWaveform_DBF(:,3);
        test_wf_freq = fftshift(fft(3*test_wf_time));
        Fs = carrier_bw;
        freq_vec = linspace(-Fs/2,Fs/2,numel(test_wf_freq)+1);
        freq_vec = freq_vec(1:end-1);
        figure
        subplot(2,1,1)
        plot(abs(test_wf_time))
        subplot(2,1,2)
        plot(freq_vec,20*log10(abs(test_wf_freq)))
    
        samplerate = info.ResourceGrids(1).Info.SampleRate;
        nfft = info.ResourceGrids(1).Info.Nfft;
        figure;
        spectrogram(test_wf_time,ones(nfft,1),0,nfft,'centered',samplerate,'yaxis','MinThreshold',-130);
        title('Spectrogram of 5G Uplink Baseband Waveform');
    
    end
    %% Analyze for EVM
    [ref_snr, ref_snr_diff,sim_params_ref] = analyze_rx_iq_samps(rxWaveform_ADC1x,sim_params,waveform_params); 
    diag_flag = 1;
    [diag_snr, diag_snr_diff,sim_params_diag] = analyze_rx_iq_samps(rxWaveform_ADC4x,sim_params,waveform_params,diag_flag); 
    [output_snr, output_snr_diff,sim_params_out] = analyze_rx_iq_samps(rxWaveform_ADC4x,sim_params,waveform_params); 
    
end