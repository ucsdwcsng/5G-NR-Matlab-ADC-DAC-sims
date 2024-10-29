function [output_snr,ref_snr] = output_snr_calib(input_snr,seed,interp_fac)
    %% Sim params
    tic
    rng(seed)
    % To print EVM statistics, set |displayEVM| to |true|. To disable the
    % prints, set |displayEVM| to |false|. To plot EVM statistics, set
    % |plotEVM| to |true|. To disable the plots, set |plotEVM| to |false|
    displayEVM = false;
    plotEVM = false;
    % To measure EVM as defined in TS 38.101-1 (FR1) or TS 38.101-2 (FR2),
    % Annex F respectively, set |evm3GPP| to |true|. |evm3GPP| is disabled by
    % default.
    evm3GPP = false;
    single_user = true;
    enable_artificial_cfo = false;
    plot_spectrum = false;
    carrier_bw = 100; % mhz
    SNR = input_snr;
    %%
    % Create a waveform generator object,then and generate the waveform.
    % wavegen = hNRReferenceWaveformGenerator(rc);
    % [txWaveform,tmwaveinfo,resourcesinfo] = generateWaveform(wavegen,wavegen.Config.NumSubframes);
    
    %Single user uplink waveform generator
    if (single_user)
        % carrier = nrSCSCarrierConfig('NSizeGrid',100);
        carrier = nrSCSCarrierConfig('SubcarrierSpacing',30,'NSizeGrid',273);
        % carrier.SubcarrierSpacing = 30;
        bwp = nrWavegenBWPConfig('SubcarrierSpacing',carrier.SubcarrierSpacing,...
                                'NSizeBWP',273);
        cfgUL = nrULCarrierConfig( ...
            'FrequencyRange','FR1', ...
            'ChannelBandwidth',carrier_bw, ...
            'NumSubframes',1, ...
            'SCSCarriers',{carrier}, ...
            'BandwidthParts',{bwp});
        cfgUL.PUSCH{1}.NumLayers = 4;
        cfgUL.PUSCH{1}.NumAntennaPorts = 4;
        cfgUL.PUSCH{1}.RVSequence = 0;
        cfgUL.PUSCH{1}.PRBSet = 0:1:272;
    else
        carriers = {
        nrSCSCarrierConfig('SubcarrierSpacing',15,'NStartGrid',10,'NSizeGrid',100), ...
        nrSCSCarrierConfig('SubcarrierSpacing',30,'NStartGrid',0,'NSizeGrid',70)};
        bwp = {
        nrWavegenBWPConfig('BandwidthPartID',0,'SubcarrierSpacing',15,'NStartBWP',30,'NSizeBWP',80), ...
        nrWavegenBWPConfig('BandwidthPartID',1,'SubcarrierSpacing',30,'NStartBWP',0,'NSizeBWP',60)};
        pusch = {
        nrWavegenPUSCHConfig('BandwidthPartID',0,'Modulation','16QAM','SlotAllocation',0:2:9,'PRBSet',0:19,'RNTI',1,'NID',1), ...
        nrWavegenPUSCHConfig('BandwidthPartID',1,'Modulation','QPSK','RNTI',2,'NID',2,'PRBSet',50:59)};
        pucch = {nrWavegenPUCCH0Config('BandwidthPartID',1,'SlotAllocation',0:9,'PRBSet',2,'DataSourceUCI', 'PN9')};
        srs = {
        nrWavegenSRSConfig('BandwidthPartID',0,'SlotAllocation',1:2:9,'NumSRSPorts',2), ... 
        nrWavegenSRSConfig('BandwidthPartID',1,'FrequencyStart',4)};
        cfgUL = nrULCarrierConfig( ...
        'FrequencyRange','FR1', ...
        'ChannelBandwidth',carrier_bw, ...
        'NumSubframes',20, ...
        'SCSCarriers',carriers, ...
        'BandwidthParts',bwp, ...
        'PUSCH',pusch, ...
        'PUCCH',pucch, ...
        'SRS',srs);
    end
    
    [txWaveform,info] = nrWaveformGenerator(cfgUL);
    %% Impairment: In-phase and Quadrature (I/Q) Imbalance, Phase Noise and Memoryless Nonlinearity
    % This example considers the most typical impairments that distort the
    % waveform when passing through an RF transmitter or receiver: phase noise,
    % I/Q imbalance, and memoryless nonlinearity. Enable or disable impairments
    % by toggling the flags |phaseNoiseOn|, |IQImbalanceON|, and
    % |nonLinearityModelOn|.
    
    IQImbalanceON = false;
    
    %%
    % Normalize the waveform to fit the dynamic range of the nonlinearity.
    txWaveform = txWaveform/max(abs(txWaveform),[],'all');
    
    % add noise
    %%
    % The waveform consists of one frame for frequency division duplexing (FDD)
    % and two for time division duplexing (TDD). Repeat the signal twice. For
    % this example, remove the first half of the resulting waveform to avoid
    % the transient introduced by the phase noise model.
    % txWaveform = repmat(txWaveform,2,1);
    
    % txwaveform_analog = 
    
    %%
    analog_factor = 160;
    adc_oversamp_factor_trad = 4;
    adc_oversamp_factor = 16;
    
    rxWaveform_analog_afterLNA_nonoise = resample(txWaveform,analog_factor,1);
    noise_sigma = sqrt(0.5*10^(-(SNR+4)/10)); % +6 normalizes
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
    % rxWaveform=txWaveform;
    % %%
    % % Introduce I/Q imbalance. Apply a 0.2 dB amplitude imbalance and a 0.5
    % % degree phase imbalance to the waveform. You can also increase the
    % % amplitude and phase imbalances by setting |amplitudeImbalance| and
    % % |phaseImbalance| to higher values.
    % if IQImbalanceON
    %     amplitudeImbalance = 0.2;
    %     phaseImbalance = 0.5;
    %     rxWaveform_DBF = iqimbal(rxWaveform_DBF,amplitudeImbalance,phaseImbalance);
    % end
    % 
    
    %% Measurements
    % The function, hNRPUSCHEVM, performs these steps to decode and analyze the
    % waveform:
    %
    % * Coarse frequency offset estimation and correction
    % * Integer frequency offset estimation and correction
    % * I/Q imbalance estimation and correction
    % * Synchronization using the demodulation reference signal (DM-RS) over
    % one frame for FDD (two frames for TDD)
    % * OFDM demodulation of the received waveform
    % * Fine frequency offset estimation and correction
    % * Channel estimation
    % * Equalization
    % * Common Phase error (CPE) estimation and compensation
    % * PUSCH EVM computation (enable the switch |evm3GPP| to process according
    % to the EVM measurement requirements specified in TS 38.101-1(FR1) or TS
    % 38.101-2 (FR2), Annex F.
    %
    % The example measures and outputs various EVM related statistics per
    % symbol, per slot, and per frame peak EVM and RMS EVM. The example
    % displays the EVM for each slot and frame, and displays the overall EVM
    % averaged over the entire input waveform. The example produces these
    % plots: EVM versus per OFDM symbol, slot, subcarrier, overall EVM and
    % in-band emissions. Each EVM related plot displays the peak versus RMS
    % EVM.
    
    % Compute and display EVM measurements
    cfg = struct();
    cfg.Evm3GPP = evm3GPP;
    cfg.PlotEVM = plotEVM;
    cfg.DisplayEVM = displayEVM;
    cfg.IQImbalance = IQImbalanceON;
    cfg.CorrectCoarseFO=true;
    cfg.CorrectFineFO=true;
    % cfgUL=wavegen.Config;
    [evmInfo1x,eqSym1x,refSym1x,hest1x,pilot_indices1x] = hNRPUSCHEVM(cfgUL,rxWaveform_ADC1x,cfg);
    [evmInfo4x,eqSym4x,refSym4x,hest4x,pilot_indices4x] = hNRPUSCHEVM(cfgUL,rxWaveform_ADC4x,cfg);
    

    
    %% Plotting
    hmat = hest4x{1};
    plot_chan = true;
    if(plot_chan)
        figure(10)
        subplot(2,1,1)
        hold on
        plot(20*log10(abs(hmat(:,1,1,1))))
        plot(20*log10(abs(hmat(:,1,1,2))))
        plot(20*log10(abs(hmat(:,1,1,3))))
        plot(20*log10(abs(hmat(:,1,1,4))))
        ylim([-40,40]);
        ylabel('Magnitude(dB)')
        legend(["H11","H12","H13","H14"])
        title('Channel magnitude(dB) and phase(degree) for Rx1')
        
        subplot(2,1,2)
        plot(180/pi*(angle(hmat(:,1,1,1))))
        hold on
        plot(180/pi*(angle(hmat(:,1,1,2))))
        plot(180/pi*(angle(hmat(:,1,1,3))))
        plot(180/pi*(angle(hmat(:,1,1,4))))
        ylim([-180,180]);
        ylabel('Angle')
        legend(["H11","H12","H13","H14"])
    end
    
    evm_abs = 10*log10(1./abs(evmInfo4x.OverallEVM.EV).^2);
    evm_snr_per_stream = mean(evm_abs,1);
    output_snr = mean(evm_snr_per_stream);

    evm_abs1x = 10*log10(1./abs(evmInfo1x.OverallEVM.EV).^2);
    evm_snr_per_stream1x = mean(evm_abs1x,1);
    ref_snr = mean(evm_snr_per_stream1x);
    

    % [evmInfo1x,eqSym1x,refSym1x,hest1x,pilot_indices1x] = hNRPUSCHEVM(cfgUL,rxWaveform_ADC1x,cfg);

    toc
    % figure()
    % subplot(2,1,1)
    % hold on
    % plot(20*log10(abs(hmat(1,:,2,1))))
    % plot(20*log10(abs(hmat(1,:,2,2))))
    % plot(20*log10(abs(hmat(1,:,2,3))))
    % plot(20*log10(abs(hmat(1,:,2,4))))
    % ylim([-40,40]);
    % ylabel('Magnitude(dB)')
    % legend(["H21","H22","H23","H24"])
    % title('Channel magnitude(dB) and phase(degree) for Rx2')
    % 
    % subplot(2,1,2)
    % plot(180/pi*(angle(hmat(1,:,1,1))))
    % hold on
    % plot(180/pi*(angle(hmat(1,:,1,2))))
    % plot(180/pi*(angle(hmat(1,:,1,3))))
    % plot(180/pi*(angle(hmat(1,:,1,4))))
    % ylim([-180,180]);
    % ylabel('Angle')
    % legend(["H21","H22","H23","H24"])
    % 
    % figure()
    % subplot(2,1,1)
    % hold on
    % plot(20*log10(abs(hmat(1,:,3,1))))
    % plot(20*log10(abs(hmat(1,:,3,2))))
    % plot(20*log10(abs(hmat(1,:,3,3))))
    % plot(20*log10(abs(hmat(1,:,3,4))))
    % ylim([-40,40]);
    % ylabel('Magnitude(dB)')
    % legend(["H31","H32","H33","H34"])
    % title('Channel magnitude(dB) and phase(degree) for Rx3')
    % 
    % subplot(2,1,2)
    % plot(180/pi*(angle(hmat(1,:,3,1))))
    % hold on
    % plot(180/pi*(angle(hmat(1,:,3,2))))
    % plot(180/pi*(angle(hmat(1,:,3,3))))
    % plot(180/pi*(angle(hmat(1,:,3,4))))
    % ylim([-180,180]);
    % ylabel('Angle')
    % legend(["H31","H32","H33","H34"])
    % 
    % figure()
    % subplot(2,1,1)
    % hold on
    % plot(20*log10(abs(hmat(1,:,4,1))))
    % plot(20*log10(abs(hmat(1,:,4,2))))
    % plot(20*log10(abs(hmat(1,:,4,3))))
    % plot(20*log10(abs(hmat(1,:,4,4))))
    % ylim([-40,40]);
    % ylabel('Magnitude(dB)')
    % legend(["H41","H42","H43","H44"])
    % title('Channel magnitude(dB) and phase(degree) for Rx4')
    % 
    % subplot(2,1,2)
    % plot(180/pi*(angle(hmat(1,:,4,1))))
    % hold on
    % plot(180/pi*(angle(hmat(1,:,4,2))))
    % plot(180/pi*(angle(hmat(1,:,4,3))))
    % plot(180/pi*(angle(hmat(1,:,4,4))))
    % ylim([-180,180]);
    % ylabel('Angle')
    % legend(["H41","H42","H43","H44"])
end