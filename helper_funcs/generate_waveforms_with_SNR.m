function waveforms = generate_waveforms_with_SNR(SNR)
    single_user = true;
    carrier_bw = 100; % mhz
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
            'NumSubframes',10, ...
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
    
    
    %%
    % Normalize the waveform to fit the dynamic range of the nonlinearity.
    txWaveform = txWaveform/max(abs(txWaveform),[],'all');
    
    % add noise
    %%
    % The waveform consists of one frame for frequency division duplexing (FDD)
    % and two for time division duplexing (TDD). Repeat the signal twice. For
    % this example, remove the first half of the resulting waveform to avoid
    % the transient introduced by the phase noise model.
    rxWaveform = repmat(txWaveform,2,1);
    
    %%
    % rxWaveform = zeros(size(txWaveform));
    % % rxWaveform(:,1) = txWaveform(:,1)+0.1*exp(1j*pi/2)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
    % % rxWaveform(:,2) = txWaveform(:,2)+0.1*exp(1j*pi/18)*txWaveform(:,1)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4);
    % % rxWaveform(:,3) = txWaveform(:,3)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,1)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
    % % rxWaveform(:,4) = txWaveform(:,4)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,1); 
    % 
    % rxWaveform(:,1) = txWaveform(:,1); 
    % rxWaveform(:,2) = txWaveform(:,2);
    % rxWaveform(:,3) = txWaveform(:,3); 
    % rxWaveform(:,4) = txWaveform(:,4); 
    
    
    noise_sigma = sqrt(0.5*10^(-(SNR+11)/10)); % +11 normalizes
    rx_waveform_size = size(rxWaveform);
    complex_noise =  noise_sigma*(randn(rx_waveform_size)+1j*randn(rx_waveform_size));
    
    rxWaveform = rxWaveform + complex_noise;

    waveforms = rxWaveform;
end