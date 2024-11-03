clear all;
%% Sim params
tic
rng(1)
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
SNR = 10;
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
        'NumSubframes',2, ...
        'SCSCarriers',{carrier}, ...
        'BandwidthParts',{bwp});
    cfgUL.PUSCH{1}.NumLayers = 4;
    cfgUL.PUSCH{1}.NumAntennaPorts = 4;
    cfgUL.PUSCH{1}.RVSequence = 0;
    cfgUL.PUSCH{1}.PRBSet = 0:1:272;
    % cfgUL.PUSCH{1}.Modulation = "16QAM";
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

%%
rxWaveform = zeros(size(txWaveform));
% rxWaveform(:,1) = txWaveform(:,1)+0.1*exp(1j*pi/2)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
% rxWaveform(:,2) = txWaveform(:,2)+0.1*exp(1j*pi/18)*txWaveform(:,1)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,4);
% rxWaveform(:,3) = txWaveform(:,3)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,1)+0.001*exp(1j*pi/36)*txWaveform(:,4); 
% rxWaveform(:,4) = txWaveform(:,4)+0.1*exp(1j*pi/18)*txWaveform(:,2)+0.01*exp(1j*pi/9)*txWaveform(:,3)+0.001*exp(1j*pi/36)*txWaveform(:,1); 

rxWaveform(:,1) = txWaveform(:,1); 
rxWaveform(:,2) = txWaveform(:,2);
rxWaveform(:,3) = txWaveform(:,3); 
rxWaveform(:,4) = txWaveform(:,4); 


noise_sigma = sqrt(0.5*10^(-(SNR+11)/10)); % +11 normalizes
rx_waveform_size = size(rxWaveform);
complex_noise =  noise_sigma*(randn(rx_waveform_size)+1j*randn(rx_waveform_size));

rxWaveform = rxWaveform + complex_noise;

if(enable_artificial_cfo)
    f_offset = 1e-3; % can go as high as 1e-3
    rxWaveform = rxWaveform.*transpose(exp(-1i*2*pi*f_offset*(0:length(rxWaveform)-1)));
end

if(plot_spectrum)
    test_wf_time = rxWaveform(:,3);
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
%%
% Introduce I/Q imbalance. Apply a 0.2 dB amplitude imbalance and a 0.5
% degree phase imbalance to the waveform. You can also increase the
% amplitude and phase imbalances by setting |amplitudeImbalance| and
% |phaseImbalance| to higher values.
if IQImbalanceON
    amplitudeImbalance = 0.2;
    phaseImbalance = 0.5;
    rxWaveform = iqimbal(rxWaveform,amplitudeImbalance,phaseImbalance);
end


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
[evmInfo,eqSym,refSym,hest,pilot_indices] = hNRPUSCHEVM(cfgUL,rxWaveform,cfg);

%% Plotting
hmat = hest{1};

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

evm_abs = 10*log10(1./abs(evmInfo.OverallEVM.EV).^2);
evm_snr_per_stream = mean(evm_abs,1)

toc

%%
% const_ref = complex(refSym{1});
% const_equalized = complex(eqSym{1});
% 
% figure
% scatter(real(const_equalized),imag(const_equalized))
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
