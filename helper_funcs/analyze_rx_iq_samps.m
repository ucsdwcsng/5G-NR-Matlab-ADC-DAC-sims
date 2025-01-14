function [mean_SNR,SNR_diff,sim_params_afterrx] = analyze_rx_iq_samps(rxWaveform,sim_params,waveform_params,diag_flag)
    
    if(nargin==3)
        diag_flag = 0;
    end

    simLocal = sim_params;
    carrier = simLocal.Carrier;
    pusch = simLocal.PUSCH;
    puschextra = simLocal.PUSCHExtension;
    correctCoarseFO_beforetimesych = false;
    correctCoarseFO_aftertimesych = true;
    sampleRate = waveform_params.SampleRate;
    % sampleRate = 20.48e6;
        % Create UL-SCH encoder System object to perform transport channel encoding
    encodeULSCH = nrULSCH;
    encodeULSCH.MultipleHARQProcesses = true;
    decodeULSCH = nrULSCHDecoder;
    decodeULSCH.MultipleHARQProcesses = true;
    decodeULSCH.TargetCodeRate = simLocal.PUSCHExtension.TargetCodeRate;
    decodeULSCH.LDPCDecodingAlgorithm = simLocal.PUSCHExtension.LDPCDecodingAlgorithm;
    decodeULSCH.MaximumLDPCIterationCount = simLocal.PUSCHExtension.MaximumLDPCIterationCount;
    pathFilters = [];

    % Create PUSCH object configured for the non-codebook transmission
    % scheme, used for receiver operations that are performed with respect
    % to the PUSCH layers
    puschNonCodebook = pusch;
    puschNonCodebook.TransmissionScheme = 'nonCodebook';

     % Specify the fixed order in which we cycle through the HARQ process IDs
    harqSequence = 0:puschextra.NHARQProcesses-1;
        
    % Set up redundancy version (RV) sequence for all HARQ processes
    if simLocal.PUSCHExtension.EnableHARQ
        % From PUSCH demodulation requirements in RAN WG4 meeting #88bis (R4-1814062)
        rvSeq = [0 2 3 1];
    else
        % HARQ disabled - single transmission with RV=0, no retransmissions
        rvSeq = 0;
    end
    
    % Initialize the state of all HARQ processes
    harqEntity = HARQEntity(harqSequence,rvSeq);

    % Calculate the transport block size for the transmission in the slot
    [puschIndices,puschIndicesInfo] = nrPUSCHIndices(carrier,pusch);
    MRB = numel(puschIndicesInfo.PRBSet);
    trBlkSize = nrTBS(pusch.Modulation,pusch.NumLayers,MRB,puschIndicesInfo.NREPerPRB,puschextra.TargetCodeRate,puschextra.XOverhead);
    
    decodeULSCHLocal = decodeULSCH;  % Copy of the decoder handle to help PCT classification of variable
    decodeULSCHLocal.reset();        % Reset decoder at the start of each SNR point
    % HARQ processing
    % If new data for current process then create a new UL-SCH transport block
    if harqEntity.NewData 
        trBlk = randi([0 1],trBlkSize,1);
        setTransportBlock(encodeULSCH,trBlk,harqEntity.HARQProcessID);
        % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
        if harqEntity.SequenceTimeout
            resetSoftBuffer(decodeULSCHLocal,harqEntity.HARQProcessID);
        end
    end

    % Encode the UL-SCH transport block
    codedTrBlock = encodeULSCH(pusch.Modulation,pusch.NumLayers, ...
        puschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

    % Create resource grid for a slot
    puschGrid = nrResourceGrid(carrier,simLocal.NTxAnts);
       
    % PUSCH modulation, including codebook based MIMO precoding if TxScheme = 'codebook'
    puschSymbols = nrPUSCH(carrier,pusch,codedTrBlock);

    % Implementation-specific PUSCH MIMO precoding and mapping. This 
    % MIMO precoding step is in addition to any codebook based 
    % MIMO precoding done during PUSCH modulation above
    if (strcmpi(pusch.TransmissionScheme,'codebook'))
        % Codebook based MIMO precoding, F precodes between PUSCH
        % transmit antenna ports and transmit antennas
        F = eye(pusch.NumAntennaPorts,simLocal.NTxAnts);
    else
        % Non-codebook based MIMO precoding, F precodes between PUSCH 
        % layers and transmit antennas
        F = eye(pusch.NumLayers,simLocal.NTxAnts);
    end
    [~,puschAntIndices] = nrExtractResources(puschIndices,puschGrid);
    puschGrid(puschAntIndices) = puschSymbols * F;
    
    % Implementation-specific PUSCH DM-RS MIMO precoding and mapping.
    % The first DM-RS creation includes codebook based MIMO precoding if applicable
    dmrsSymbols = nrPUSCHDMRS(carrier,pusch);
    dmrsIndices = nrPUSCHDMRSIndices(carrier,pusch);
    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
        puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
    end
    
    % CFO correction (coarse)
    % 
    % foffsetEstCoarse = 0;
    % integerCfoHz = 0;
    % 
    % [foffsetEstCoarse,offset] = hNRFrequencyOffset('coarseFOwoffset',carrier,rxWaveform,sampleRate);
    % 
    % if correctCoarseFO_beforetimesych
    %     rxWaveform= hNRFrequencyOffset('FOCorrect',rxWaveform,sampleRate,foffsetEstCoarse);
    % 
    %     % integerCfoHz = hNRFrequencyOffset('integerFO',carrier,rxWaveform,refGrid,sampleRate);
    %     % rxWaveform = hNRFrequencyOffset('FOCorrect',rxWaveform,sampleRate,integerCfoHz);
    % end

    % Practical synchronization. Correlate the received waveform 
    % with the PUSCH DM-RS to give timing offset estimate 't' and
    % correlation magnitude 'mag'. The function
    % hSkipWeakTimingOffset is used to update the receiver timing
    % offset. If the correlation peak in 'mag' is weak, the current
    % timing estimate 't' is ignored and the previous estimate
    % 'offset' is used
    [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
    [corr_val,offset] = max(mag(1:simLocal.txslotsize_samples+simLocal.num_zeros_to_prepend));
    % offset = 43380;
    % offset = 26686;
    % offset = offset;
    
    % offset = hSkipWeakTimingOffset(0,t,mag());
    % Display a warning if the estimated timing offset exceeds the
    % maximum channel delay
    % if offset > maxChDelay
    %     warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
    %         ' This will result in a decoding failure. This may be caused by low SNR,' ...
    %         ' or not enough DM-RS symbols to synchronize successfully.'],offset,maxChDelay);
    % end
    % + 100 leeway
    nSlots = floor((size(rxWaveform,1)-offset+simLocal.num_zeros_to_prepend)/(simLocal.txslotsize_samples+simLocal.num_zeros_to_prepend));
    slotsize_samps = simLocal.txslotsize_samples;
    
    nslot_index = 1;
    
    % nSlots = 3;
    slot_offset = 0;

    for nslot = 0:1:nSlots-1
        % nSlots 
        rxWaveform_to_analyze = rxWaveform(offset+(nslot+slot_offset)*(slotsize_samps+simLocal.num_zeros_to_prepend):...
            offset+(nslot+slot_offset)*(slotsize_samps+simLocal.num_zeros_to_prepend)+slotsize_samps,:);
        
        foffsetEstCoarse = 0;
        integerCfoHz = 0;
        if correctCoarseFO_aftertimesych
            %   FOFFSET =
            %   hNRFrequencyOffset('COARSEFO',CARRIER,WAVEFORM,SAMPLERATE,TOFFSET)
            % foffsetEstCoarse = hNRFrequencyOffset('coarseFO',carrier,rxWaveform,sampleRate,offset)
            foffsetEstCoarse = hNRFrequencyOffset('coarseFO',carrier,rxWaveform_to_analyze,sampleRate);
            rxWaveform_to_analyze = hNRFrequencyOffset('FOCorrect',rxWaveform_to_analyze,sampleRate,foffsetEstCoarse);
    
            % integerCfoHz = hNRFrequencyOffset('integerFO',carrier,rxWaveform,refGrid,sampleRate);
            % rxWaveform = hNRFrequencyOffset('FOCorrect',rxWaveform,sampleRate,integerCfoHz);
        end

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform_to_analyze);
        [K,L,R] = size(rxGrid);
        if (L < carrier.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
        end
            
        % if (simLocal.PerfectChannelEstimator)
        %     % For perfect channel estimate, use the OFDM channel response
        %     % obtained from the channel
        %     estChannelGrid = ofdmResponse;
        % 
        %     % For perfect noise estimate, use the noise power per RE
        %     noiseEst = nPowerPerRE;
        % 
        %     % Apply MIMO deprecoding to estChannelGrid to give an estimate
        %     % per transmission layer
        %     K = size(estChannelGrid,1);
        %     estChannelGrid = reshape(estChannelGrid,K*carrier.SymbolsPerSlot*simLocal.NRxAnts,simLocal.NTxAnts);
        %     estChannelGrid = estChannelGrid * F.';
        %     if (strcmpi(pusch.TransmissionScheme,'codebook'))
        %         W = nrPUSCHCodebook(pusch.NumLayers,pusch.NumAntennaPorts,pusch.TPMI,pusch.TransformPrecoding);
        %         estChannelGrid = estChannelGrid * W.';
        %     end
        %     estChannelGrid = reshape(estChannelGrid,K,carrier.SymbolsPerSlot,simLocal.NRxAnts,[]);
        % else
            
        % Practical channel estimation between the received grid and
        % each transmission layer, using the PUSCH DM-RS for each layer
        % which are created by specifying the non-codebook transmission
        % scheme
        dmrsLayerSymbols = nrPUSCHDMRS(carrier,puschNonCodebook);
        dmrsLayerIndices = nrPUSCHDMRSIndices(carrier,puschNonCodebook);
        % [estChannelGrid,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsLayerIndices,dmrsLayerSymbols,'CDMLengths',pusch.DMRS.CDMLengths);
        [estChannelGrid,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsLayerIndices,dmrsLayerSymbols);
        % end
            
        % Get PUSCH resource elements from the received grid
        [puschRx,puschHest] = nrExtractResources(puschIndices,rxGrid,estChannelGrid);
        
        channel_vis = false;
        if(channel_vis)
            hmat = puschHest;
            sensitivity_offset = 270;
            figure(10)
            subplot(2,1,1)
            hold on
            plot(20*log10(abs(hmat(:,1,1)))+sensitivity_offset)
            plot(20*log10(abs(hmat(:,1,2)))+sensitivity_offset)
            plot(20*log10(abs(hmat(:,2,1)))+sensitivity_offset)
            plot(20*log10(abs(hmat(:,2,2)))+sensitivity_offset)
            % plot(20*log10(abs(hmat(:,1,3))))
            % plot(20*log10(abs(hmat(:,1,4))))
            % plot(20*log10(abs(hmat(:,1)))+230)
            % plot(20*log10(abs(hmat(:,1))))
            % plot(20*log10(abs(hmat(:,1))))
            % plot(20*log10(abs(hmat(:,1))))
            % ylim([-40,40]);
            ylabel('Magnitude(dB)')
            legend(["H11","H12","H21","H22"])
            title('Channel magnitude(dB) and phase(degree) for Rx1')
            
            subplot(2,1,2)
            % plot(180/pi*(unwrap(angle(hmat(:,1)))))
            plot(180/pi*(unwrap(angle(hmat(:,1,1)))))
            hold on
            plot(180/pi*(unwrap(angle(hmat(:,2,2)))))
            legend(["H11","H22"])
            
            ylim([-220,220]);
            ylabel('Angle')

            differential_phase = mean(180/pi*(unwrap(angle(hmat(:,1,1)))) - 180/pi*(unwrap(angle(hmat(:,2,2)))))

            % legend(["H11","H12","H13","H14"])
        end
    
        % Equalization
        % [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
        puschHest_diag = zeros(size(puschHest));
        puschHest_diag(:,1,1) =  puschHest(:,1,1);
        puschHest_diag(:,2,2) =  puschHest(:,2,2);
        puschHest_diag(:,3,3) =  puschHest(:,3,3);
        puschHest_diag(:,4,4) =  puschHest(:,4,4);
        
        if(diag_flag==1)
            [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest_diag,0);
        else
            [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest,0);
        end
        
        
        % Display EVM per layer, per slot and per RB. Reference symbols for
        % each layer are created by specifying the non-codebook
        % transmission scheme
        
        % refSymbols = nrPUSCH(carrier,puschNonCodebook,codedTrBlock);
        % meanSNR = plotLayerEVM(nSlots,nslot,puschNonCodebook,size(puschGrid),puschIndices,simLocal.refSymbols{nslot_index},puschEq);
        [mean_SNR,SNR_diff] = calcEVM(puschNonCodebook,size(puschGrid),puschIndices,simLocal.refSymbols{nslot_index},puschEq);
        % meanSNR = plotLayerEVM(nSlots,nslot,puschNonCodebook,size(puschGrid),puschIndices,simLocal.refSymbols{nslot_index},rxSymbols)
        % mean_SNR_per_slot = meanSNR;

        fprintf('\n(%3.2f%%) NSlot=%d, Processed\n',100*(nslot+1)/nSlots,nslot);
        %% if using multi slot witthin the protocol
        %% NOT RECOMMENDED FOr HWR
        if(simLocal.num_slots_to_gen>1)
            nslot_index = mod(nslot_index+1,simLocal.num_slots_to_gen);
        end

    end
    

    % fprintf('\nThroughput(Mbps) for %d frame(s) = %.4f\n',simLocal.NFrames,1e-6*simThroughput(snrIdx)/(simLocal.NFrames*10e-3));
    % fprintf('Throughput(%%) for %d frame(s) = %.4f\n',simLocal.NFrames,simThroughput(snrIdx)*100/maxThroughput(snrIdx));
    sim_params_afterrx = simLocal;
end

function mean_SNR = plotLayerEVM(NSlots,nslot,pusch,siz,puschIndices,puschSymbols,puschEq)
% Plot EVM information

    persistent slotEVM;
    persistent rbEVM
    persistent evmPerSlot;
    
    if (nslot==0)
        slotEVM = comm.EVM;
        rbEVM = comm.EVM;
        evmPerSlot = NaN(NSlots,pusch.NumLayers);
        
    end
    figure(1);
    evmPerSlot(nslot+1,:) = slotEVM(puschSymbols,puschEq);
    subplot(2,1,1);
    plot(0:(NSlots-1),evmPerSlot,'o-');
    xlabel('Slot number');
    ylabel('EVM (%)');
    legend("layer " + (1:pusch.NumLayers),'Location','EastOutside');
    title('EVM per layer per slot');
    hold on

    subplot(2,1,2);
    [k,~,p] = ind2sub(siz,puschIndices);
    rbsubs = floor((k-1) / 12);
    NRB = siz(1) / 12;
    evmPerRB = NaN(NRB,pusch.NumLayers);
    for nu = 1:pusch.NumLayers
        for rb = unique(rbsubs).'
            this = (rbsubs==rb & p==nu);
            evmPerRB(rb+1,nu) = rbEVM(puschSymbols(this),puschEq(this));
        end
    end
    plot(0:(NRB-1),evmPerRB,'x-');
    hold on
    xlabel('Resource block');
    ylabel('EVM (%)');
    legend("layer " + (1:pusch.NumLayers),'Location','EastOutside');
    title(['EVM per layer per resource block, slot #' num2str(nslot)]);
    
    drawnow;
    SNR_dB = 10*log10(1./(0.01*evmPerRB).^2)-3.7;
    mean_SNR = mean(SNR_dB,1)
    % 
    % evm_abs = 10*log10(1./abs(evmInfo.OverallEVM.EV).^2);
    % evm_snr_per_stream = mean(evm_abs,1)
end

function [mean_SNR,SNR_diff] = calcEVM(pusch,siz,puschIndices,puschSymbols,puschEq)
    rbEVM = comm.EVM;
    [k,~,p] = ind2sub(siz,puschIndices);
    rbsubs = floor((k-1) / 12);
    NRB = siz(1) / 12;
    evmPerRB = NaN(NRB,pusch.NumLayers);
    for nu = 1:pusch.NumLayers
        for rb = unique(rbsubs).'
            this = (rbsubs==rb & p==nu);
            evmPerRB(rb+1,nu) = rbEVM(puschSymbols(this),puschEq(this));
        end
    end
    
    SNR_dB = 10*log10(1./(0.01*evmPerRB).^2)-3.7;
    SNR_mean_across_subC = mean(SNR_dB,1);
    mean_SNR = mean(SNR_mean_across_subC);
    SNR_diff = max(SNR_mean_across_subC)-min(SNR_mean_across_subC);
    % 
    % evm_abs = 10*log10(1./abs(evmInfo.OverallEVM.EV).^2);
    % evm_snr_per_stream = mean(evm_abs,1)
end
