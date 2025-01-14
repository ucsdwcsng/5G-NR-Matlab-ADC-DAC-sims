function [tx_iqsamples,simParameters,waveformParams] = get_tx_iq_samples_mimo(num_slots_to_gen,SCS,grid_size,num_layers,num_tx_ant,num_rx_ant,num_zeros_to_append,num_zeros_to_prepend)
    if(nargin==6)
        num_zeros_to_append = 0;
        num_zeros_to_prepend = 1000;
    end
    % disp("Here")
    simParameters = struct();       % Clear simParameters variable to contain all key simulation parameters
    simParameters.NFrames = 1;      % Number of 10 ms frames
    
    % Set waveform type and PUSCH numerology (SCS and CP type)
    simParameters.Carrier = nrCarrierConfig;        % Carrier resource grid configuration
    simParameters.Carrier.NSizeGrid = grid_size;           % Bandwidth in number of resource blocks (52 RBs at 15 kHz SCS for 10 MHz BW)
    simParameters.Carrier.SubcarrierSpacing = SCS;   % 15, 30, 60, 120 (kHz)
    simParameters.Carrier.CyclicPrefix = 'Normal';  % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
    simParameters.Carrier.NCellID = 0;              % Cell identity
    simParameters.maxframesize_in_samples = simParameters.Carrier.NSizeGrid*12*simParameters.Carrier.SlotsPerFrame;
    % simParameters.Carrier.SlotsPerFrame = 5;              % Slots per frame
    
    % PUSCH/UL-SCH parameters
    simParameters.PUSCH = nrPUSCHConfig;      % This PUSCH definition is the basis for all PUSCH transmissions in the BLER simulation
    simParameters.PUSCHExtension = struct();  % This structure is to hold additional simulation parameters for the UL-SCH and PUSCH
    
    % Define PUSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
    simParameters.PUSCH.PRBSet =  0:simParameters.Carrier.NSizeGrid-1; % PUSCH PRB allocation
    simParameters.PUSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot]; % PUSCH symbol allocation in each slot
    simParameters.PUSCH.MappingType = 'A'; % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
    
    % Scrambling identifiers
    simParameters.PUSCH.NID = simParameters.Carrier.NCellID;
    simParameters.PUSCH.RNTI = 1;
    
    % Define the transform precoding enabling, layering and transmission scheme
    simParameters.PUSCH.TransformPrecoding = false; % Enable/disable transform precoding
    simParameters.PUSCH.NumLayers = num_layers;              % Number of PUSCH transmission layers
    simParameters.PUSCH.TransmissionScheme = 'nonCodebook'; % Transmission scheme ('nonCodebook','codebook')
    simParameters.PUSCH.TPMI = 0;                   % Precoding matrix indicator for codebook based precoding
    
    % Define codeword modulation
    simParameters.PUSCH.Modulation = 'QPSK'; % 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
    % simParameters.PUSCH.Modulation = '64QAM'; % 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
    
    % PUSCH DM-RS configuration
    simParameters.PUSCH.DMRS.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
    simParameters.PUSCH.DMRS.DMRSLength = 2;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
    simParameters.PUSCH.DMRS.DMRSAdditionalPosition = 1;  % Additional DM-RS symbol positions (max range 0...3)
    simParameters.PUSCH.DMRS.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
    simParameters.PUSCH.DMRS.NumCDMGroupsWithoutData = 1; % Number of CDM groups without data
    simParameters.PUSCH.DMRS.NIDNSCID = 10;                % Scrambling identity (0...65535)
    simParameters.PUSCH.DMRS.NSCID = 1;                   % Scrambling initialization (0,1)
    simParameters.PUSCH.DMRS.NRSID = 1;                   % Scrambling ID for low-PAPR sequences (0...1007)
    simParameters.PUSCH.DMRS.GroupHopping = 0;            % Group hopping (0,1)
    simParameters.PUSCH.DMRS.SequenceHopping = 0;         % Sequence hopping (0,1)
    % 
    % simParameters.PUSCH.DMRS.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
    % simParameters.PUSCH.DMRS.DMRSLength = 1;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
    % simParameters.PUSCH.DMRS.DMRSAdditionalPosition = 1;  % Additional DM-RS symbol positions (max range 0...3)
    % simParameters.PUSCH.DMRS.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
    % simParameters.PUSCH.DMRS.NumCDMGroupsWithoutData = 1; % Number of CDM groups without data
    % simParameters.PUSCH.DMRS.NIDNSCID = 10;                % Scrambling identity (0...65535)
    % simParameters.PUSCH.DMRS.NSCID = 1;                   % Scrambling initialization (0,1)
    % simParameters.PUSCH.DMRS.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
    % simParameters.PUSCH.DMRS.GroupHopping = 0;            % Group hopping (0,1)
    % simParameters.PUSCH.DMRS.SequenceHopping = 0;         % Sequence hopping (0,1)
    % Additional simulation and UL-SCH related parameters
    %
    % Target code rate
    simParameters.PUSCHExtension.TargetCodeRate = 193 / 1024; % Code rate used to calculate transport block size
    %
    % HARQ process and rate matching/TBS parameters
    simParameters.PUSCHExtension.XOverhead = 0;       % Set PUSCH rate matching overhead for TBS (Xoh)
    simParameters.PUSCHExtension.NHARQProcesses = 16; % Number of parallel HARQ processes to use
    simParameters.PUSCHExtension.EnableHARQ = false;   % Enable retransmissions for each process, using RV sequence [0,2,3,1]
    
    % LDPC decoder parameters
    % Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
    simParameters.PUSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
    simParameters.PUSCHExtension.MaximumLDPCIterationCount = 6;
    
    % Define the overall transmission antenna geometry at end-points
    % If using a CDL propagation channel then the integer number of antenna elements is
    % turned into an antenna panel configured when the channel model object is created
    simParameters.NTxAnts = num_tx_ant; % Number of transmit antennas
    simParameters.NRxAnts = num_rx_ant; % Number of receive antennas
    simParameters.PUSCH.NumAntennaPorts = num_layers;        % Number of antenna ports for codebook based precoding
    % % Define the general CDL/TDL propagation channel parameters
    % simParameters.DelayProfile = 'TDL-A'; % Use TDL-A model (Indoor hotspot model)
    % simParameters.DelaySpread = 30e-9;
    % simParameters.MaximumDopplerShift = 10;
    % 
    % Cross-check the PUSCH layering against the channel geometry 
    validateNumLayers(simParameters);
    
    %%
    % The simulation relies on various pieces of information about the baseband 
    % waveform, such as sample rate.
    
    waveformInfo = nrOFDMInfo(simParameters.Carrier); % Get information about the baseband waveform after OFDM modulation step
    
    % Set up redundancy version (RV) sequence for all HARQ processes
    if simParameters.PUSCHExtension.EnableHARQ
        % From PUSCH demodulation requirements in RAN WG4 meeting #88bis (R4-1814062)
        rvSeq = [0 2 3 1];
    else
        % HARQ disabled - single transmission with RV=0, no retransmissions
        rvSeq = 0;
    end
    
    % Create UL-SCH encoder System object to perform transport channel encoding
    encodeULSCH = nrULSCH;
    encodeULSCH.MultipleHARQProcesses = true;
    encodeULSCH.TargetCodeRate = simParameters.PUSCHExtension.TargetCodeRate;

    decodeULSCH = nrULSCHDecoder;
    decodeULSCH.MultipleHARQProcesses = true;
    decodeULSCH.TargetCodeRate = simParameters.PUSCHExtension.TargetCodeRate;
    decodeULSCH.LDPCDecodingAlgorithm = simParameters.PUSCHExtension.LDPCDecodingAlgorithm;
    decodeULSCH.MaximumLDPCIterationCount = simParameters.PUSCHExtension.MaximumLDPCIterationCount;
    % Take copies of channel-level parameters to simplify subsequent parameter referencing 
    simLocal = simParameters;
    carrier = simLocal.Carrier;
    pusch = simLocal.PUSCH;
    puschextra = simLocal.PUSCHExtension;
    
    pathFilters = [];

    % Create PUSCH object configured for the non-codebook transmission
    % scheme, used for receiver operations that are performed with respect
    % to the PUSCH layers
    puschNonCodebook = pusch;
    puschNonCodebook.TransmissionScheme = 'nonCodebook';

     % Specify the fixed order in which we cycle through the HARQ process IDs
    harqSequence = 0:puschextra.NHARQProcesses-1;

    % Initialize the state of all HARQ processes
    harqEntity = HARQEntity(harqSequence,rvSeq);

    
    % Total number of slots in the simulation period
    % NSlots = simLocal.NFrames * carrier.SlotsPerFrame;
    NSlots = num_slots_to_gen;
    simParameters.num_slots_to_gen = num_slots_to_gen;
    
    tx_iqsamples = [];

    for nslot = 0:NSlots-1
        % Calculate the transport block size for the transmission in the slot
        [puschIndices,puschIndicesInfo] = nrPUSCHIndices(carrier,pusch);
        MRB = numel(puschIndicesInfo.PRBSet);
        trBlkSize = nrTBS(pusch.Modulation,pusch.NumLayers,MRB,puschIndicesInfo.NREPerPRB,puschextra.TargetCodeRate,puschextra.XOverhead);
    
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
        
        % OFDM modulation
        txWaveform = nrOFDMModulate(carrier,puschGrid);
        txWaveform_withzeros = padarray(txWaveform,[num_zeros_to_prepend,0],0,'pre');
        simParameters.num_zeros_to_prepend = num_zeros_to_prepend;
        % txWaveform = padarray(txWaveform,[num_zeros_to_append,0],0,'post');
        simParameters.txslotsize_samples = size(txWaveform,1);
        tx_iqsamples = [tx_iqsamples txWaveform_withzeros.'];
        simParameters.refSymbols{nslot+1} = nrPUSCH(carrier,puschNonCodebook,codedTrBlock);
    end

    % tx_iqsamples = [tx_iqsamples zeros([1, 10000])];
    waveformParams = waveformInfo;

    tx_iqsamples = 0.8*(tx_iqsamples./max(abs(tx_iqsamples),[],2));
    
end

