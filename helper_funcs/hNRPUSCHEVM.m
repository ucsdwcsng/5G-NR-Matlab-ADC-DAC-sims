function [evmInfo,eqSym,refSym,hest,refGrid_indices] = hNRPUSCHEVM(varargin) %added hest,ref_indices as output to extract channel
%hNRPUSCHEVM EVM measurement of NR PUSCH
%   [EVMINFO,EQSYM,REFSYM] = hNRPUSCHEVM(WAVECONFIG,RXWAVEFORM)
%   Calculates the error vector magnitude (EVM) of a received waveform. If
%   cfg.EVM3GPP is true, EVM measurement is done using the 3GPP specified
%   EVM algorithm as defined in TS 38.101-2 (FR2), Annex F
%
%   EVMINFO is a numBWPs-by-1 struct array. numBWPs is the number of
%   configured bandwidth parts. EVMINFO contains EVM statistics with
%   fields: 
%     SubcarrierRMS    - RMS EVM per subcarrier
%                        (Column vector of N subcarriers)
%     SubcarrierPeak   - Peak EVM per subcarrier
%                        (Column vector of N subcarriers)
%     SymbolRMS          - RMS EVM Per symbol
%                        (Column vector of S symbols)
%     SymbolPeak       - Peak EVM Per symbol
%                        (Column vector of S symbols)
%     SlotRMS          - RMS EVM per slot
%                        (Column vector of x slots)
%     SlotPeak         - Peak EVM Per slot
%                        (Column vector of x slots)
%     EVMGrid          - Raw error vector for the selected window edge
%                        (N subcarriers-by-S symbols)
%     OverallEVM       - Structure containing EVM statistics for the
%                        overall waveform. It contains these fields:
%       EV             - Raw error vector for the overall waveform
%                        (Array of 1-by-n layers)
%       RMS            - RMS EVM for the overall waveform
%                        (Scalar)
%       Peak           - Peak EVM for the overall waveform
%                        (Scalar)
%     BandwidthPartID  - Index of the bandwidth part
%     Emissions        - structure containing in-band emissions statistics
%                        as defined in TS 38.101-1 (FR1) / TS 38.101-2
%                        (FR2), Annex F.3. It contains these fields:
%       DeltaRB        - Array of unallocated PUSCH RBs.
%                        Each RB is expressed as an offset relative to the
%                        nearest edge of the PUSCH allocation
%                        (Row vector of delta RBs)
%       Absolute       - Absolute in-band emisson in each delta RB vs slot.
%                        (Array of deltaRBs vs number of slots)
%       Relative       - Relative in-band emisson in each delta RB vs slot.
%                        (Array of deltaRBs vs number of slots)
%     AverageFO        - Estimated frequency estimate
%     AmpImbalance     - Estimated amplitude imbalance (dB)
%     PhImbalance      - Estimated phase imbalance (degrees)
%   EQSYM              - Output cell array of IQ constellations for low &
%                        high EVM window locations
%                        (2-by-numBWPs-by-number of slots)
%   REFSYM             - Output cell array of reference IQ constellations
%                        for low and high EVM window locations. 
%                        (2-by-numBWPs-by-number of slots)
%   WAVECONFIG         - Input structure or object of type
%                        'nrULCarrierConfig', containing carrier and cell
%                        related parameters
%   RXWAVEFORM         - Time domain baseband IQ samples input. Timing of
%                        the waveform is assumed to be slot-wise aligned
%                        with sample level fine-tuning performed
%                        subsequently. The length of the waveform can be an
%                        arbitrary number of slots
%
%   [EVMINFO,EQSYM,REFSYM] = hNRPUSCHEVM(...,CFG) CFG is an optional
%   structure with the fields.
%   Evm3GPP         - Enables or disables 3GPP method of EVM computation.
%                     (Default value : false)
%   PlotEVM         - Enables or disables plotting of EVM (per slot, per
%                     symbol, per subcarrier and overall EVM)
%                     (Default value : true)
%   DisplayEVM      - Enables or disables the display of
%                     EVM statistics on the command window
%                     (Default value : true)
%   InitialNSlot    - Starting slot of the input waveform
%                     (Default value : 0)
%   SampleRate      - Waveform sample rate. If absent, the default
%                     waveform sample rate specified in WAVECONFIG is
%                     used for demodulation. It is specified as either a
%                     positive scalar or [].
%   CorrectCoarseFO - Enables or disables coarse frequency offset (FO)
%                     estimation and correction. If enabled, the waveform
%                     start slot should match InitialNSlot.
%                     (Default value : false)
%   CorrectFineFO   - Enables or disables fine FO estimation and correction
%                     (Default value : false)
%   TimeSyncEnable  - Enables or disables timing estimation
%                     (Default value : true)
%   UseWholeGrid    - Enables or disables the use of a reference grid,
%                     consisting of known data and demodulation reference
%                     signals (DM-RS) IQ samples for timing estimation
%                     purposes. Enable 'UseWholeGrid' when there are only
%                     few DM-RS symbols in the received waveform or when
%                     the timing synchronization is done over a short part
%                     of the frame, for example, just a slot
%                     (Default value : false)
%   IQImbalance     - Enables or disables I/Q imbalance estimation and
%                     correction
%                     (Default value : false)

% Copyright 2020-2022 The MathWorks, Inc.

    narginchk(2,3);

    waveConfig = varargin{1};
    rxWaveform = varargin{2};
    cfg = [];
    if nargin == 3
        cfg = varargin{3};
    end

    % Validate cfg
    validateInputs(cfg);

    if ~isfield(cfg,'Evm3GPP')
        evm3GPP = false;
    else
        evm3GPP = cfg.Evm3GPP;
    end
    if ~isfield(cfg,'PlotEVM')
        plotEVM = true;
    else
        plotEVM = cfg.PlotEVM;
    end
    if ~isfield(cfg,'DisplayEVM')
        displayEVM = true;
    else
        displayEVM = cfg.DisplayEVM;
    end
    if ~isfield(cfg,'InitialNSlot')
        initialNSlot = 0;
    else
        initialNSlot = cfg.InitialNSlot;
    end
    if ~isfield(cfg,'CorrectCoarseFO')
        correctCoarseFO = false;
    else
        correctCoarseFO = cfg.CorrectCoarseFO;
    end
    if ~isfield(cfg,'CorrectFineFO')
        correctFineFO = false;
    else
        correctFineFO = cfg.CorrectFineFO;
    end
    if ~isfield(cfg,'TimeSyncEnable')
        timeSyncEnable = true;
    else
        timeSyncEnable = cfg.TimeSyncEnable;
    end
    if ~isfield(cfg,'UseWholeGrid')
        useWholeGrid = false;
    else
        useWholeGrid = cfg.UseWholeGrid;
    end
    if ~isfield(cfg,'IQImbalance')
        iqImbalance = false;
    else
        iqImbalance = true;
    end

    % Get per-slot PUSCH resources (waveformResources) used as reference for EVM calculation
    % waveConfig can be an object of type 'nrULCarrierConfig' or a struct.

    if isa(waveConfig,'nrULCarrierConfig')
        [~,winfo] = nrWaveformGenerator(waveConfig);
    else
        % Convert struct to an object of type 'nrULCarrierConfig'
        waveConfig = mapStruct2ObjUplink(waveConfig);
        [~,winfo] = nrWaveformGenerator(waveConfig);
    end
    waveformResources = winfo.WaveformResources;

    % Obtain number of bandwidth parts (BWP)
    numBWPs = length(waveConfig.BandwidthParts);

    eqSym = cell(2,numBWPs);                  % Equalized symbols for constellation plot, for each low/high EVM window location and BWP
    refSym = cell(2,numBWPs);                 % Reference symbols for constellation plot, for each low/high EVM window location and BWP
    evmInfo = [];                             % EVM statistics

    % Store the received waveform for reuse in each BWP
    % Loop over each BWP
    rxWaveformOrig = rxWaveform;
    for bwpIdx = 1:numBWPs
        % Check validity of each BWP configuration
        invalidBWPConfigFlag = false;

        % Get a copy of the carrier config associated with the BWP numerology
        carrierID = nr5g.internal.wavegen.getCarrierIDByBWPIndex(waveConfig.SCSCarriers, waveConfig.BandwidthParts, bwpIdx);

        % Extract PUSCH resources and carrier related information
        puschWaveCfg = [];
        temp = [waveConfig.PUSCH{:}];
        loc = find([temp.BandwidthPartID] == waveConfig.BandwidthParts{bwpIdx}.BandwidthPartID);
        if length(loc) > 1
            warning('Multiple PUSCHs per BWP are not supported');
            loc = loc(1);
            invalidBWPConfigFlag = true;
        end
        puschDefs = [waveConfig.PUSCH{loc}];
        gridsize = waveConfig.SCSCarriers{carrierID}.NSizeGrid;
        gridstart = waveConfig.SCSCarriers{carrierID}.NStartGrid;
        cp = waveConfig.BandwidthParts{bwpIdx}.CyclicPrefix;
        scs = waveConfig.SCSCarriers{carrierID}.SubcarrierSpacing;

        if ~isempty(loc)
            ri = waveformResources.PUSCH(loc);
            ri.PUSCH = puschDefs;
            puschWaveCfg = [puschWaveCfg ri];

            % Process only if PUSCH.Enable is true
            if ~(puschWaveCfg.PUSCH.Enable) %#ok<BDLGI,*BDSCI>
                invalidBWPConfigFlag = true;
            end

            % Ensure puschWaveCfg contains the resources needed for EVM measurement
            if isempty(puschWaveCfg.Resources) || isempty(puschWaveCfg.PUSCH.PRBSet) || ...
                    isempty(puschWaveCfg.PUSCH.SymbolAllocation) || puschWaveCfg.PUSCH.SymbolAllocation(2) == 0
                warning('Input configuration does not contain adequate resources to proceed with EVM measurement');
                invalidBWPConfigFlag = true;
            end
        else
            invalidBWPConfigFlag = true;
        end

        % Skip this BWP due to unexpected configuration
        if invalidBWPConfigFlag
            evmInfo(bwpIdx).SubcarrierRMS = [];
            evmInfo(bwpIdx).SubcarrierPeak = [];
            evmInfo(bwpIdx).SymbolRMS = [];
            evmInfo(bwpIdx).SymbolPeak = [];
            evmInfo(bwpIdx).SlotRMS = [];
            evmInfo(bwpIdx).SlotPeak = [];
            evmInfo(bwpIdx).EVMGrid = [];
            evmInfo(bwpIdx).OverallEVM.EV = [];
            evmInfo(bwpIdx).OverallEVM.RMS = [];
            evmInfo(bwpIdx).OverallEVM.Peak = [];
            evmInfo(bwpIdx).Emissions = [];
            evmInfo(bwpIdx).BandwidthPartID = [];
            evmInfo(bwpIdx).AverageFO = 0;
            evmInfo(bwpIdx).AmpImbalance = 0;
            evmInfo(bwpIdx).PhImbalance = 0;
            continue;
        end

        % Extract bandwidth part related configuration
        bwpCfg = waveConfig.BandwidthParts{bwpIdx};
        bwpStart = bwpCfg.NStartBWP-gridstart;
        bwpLen = bwpCfg.NSizeBWP;

        % Create carrier resource grid configuration for synchronization and OFDM demodulation
        carrier = nrCarrierConfig;
        carrier.NCellID = waveConfig.NCellID;
        carrier.NSizeGrid = gridsize;
        carrier.NStartGrid = gridstart;
        carrier.SubcarrierSpacing = scs;
        carrier.CyclicPrefix = cp;
        carrier.NSlot = initialNSlot;

        if isempty(puschWaveCfg.PUSCH.NID)
            puschWaveCfg.PUSCH.NID = carrier.NCellID;
        end

        % Obtain OFDM related information
        ofdmInfo = nrOFDMInfo(carrier);
        ofdmInfo.SamplesPerSubframe = sum(ofdmInfo.SymbolLengths);
        sampleRate = winfo.ResourceGrids(bwpIdx).Info.SampleRate;
        if ~isfield(cfg,'SampleRate') || isempty(cfg.SampleRate)
            rxWaveform = rxWaveformOrig;
        else
            rxWaveform = resample(rxWaveformOrig,sampleRate,cfg.SampleRate);
        end
        k0 = winfo.ResourceGrids(bwpIdx).Info.k0;
        scs = carrier.SubcarrierSpacing;
        ofdmInfo.SamplesPerSubframe = sampleRate/1000;

        % Generate a reference grid of length two frames for timing synchronization
        refGrid = hReferenceGrid(carrier,bwpCfg,puschWaveCfg,2*ofdmInfo.SlotsPerFrame,'PUSCH');
        refGrid_indices = find(refGrid~=0);
        % Shift in frequency the waveform, taking into account the 'k0' for the current BWP
        t = (0:size(rxWaveform,1)-1).'/sampleRate;
        k0Offset = k0*scs*1e3;
        rxWaveformk0Shifted = rxWaveform.*exp(-1i*2*pi*k0Offset*t);

        foffsetEstCoarse = 0;
        integerCfoHz = 0;
        if correctCoarseFO
            foffsetEstCoarse = hNRFrequencyOffset('coarseFO',carrier,rxWaveformk0Shifted,sampleRate);
            rxWaveformk0Shifted = hNRFrequencyOffset('FOCorrect',rxWaveformk0Shifted,sampleRate,foffsetEstCoarse);

            integerCfoHz = hNRFrequencyOffset('integerFO',carrier,rxWaveformk0Shifted,refGrid,sampleRate);
            rxWaveformk0Shifted = hNRFrequencyOffset('FOCorrect',rxWaveformk0Shifted,sampleRate,integerCfoHz);
        end

        % Estimate and correct the I/Q imbalance for the received waveform
        ampImbEst = 0;
        phImbEst = 0;
        if iqImbalance
            [rxWaveformk0Shifted,ampImbEst,phImbEst] = hNRIQImbalance(rxWaveformk0Shifted);
        end

        % Time synchronization of input waveform
        offset = 0;
        if timeSyncEnable
            % If enabled, use a grid of known data and DM-RS IQ samples for
            % timing estimation
            if useWholeGrid
                refGrid = winfo.ResourceGrids(bwpIdx).ResourceGridInCarrier;
                refGrid(:,1:ofdmInfo.SymbolsPerSlot*carrier.NSlot,:) = [];
            end
            if ~isempty(refGrid)
                offset = nrTimingEstimate(carrier,rxWaveformk0Shifted,refGrid,'SampleRate',sampleRate);
            end
        end
        rxWaveform = rxWaveformk0Shifted(1+offset:end,:);

        % Calculate number of slots and frames for the given sampleRate
        L = carrier.SymbolsPerSlot;
        grid = nrOFDMDemodulate(carrier,rxWaveform,'SampleRate',sampleRate);
        nSlots = floor(size(grid,2)/L);
        nFrames = floor(nSlots/(10*ofdmInfo.SlotsPerSubframe));

        % When the input waveform has more slots than
        % waveConfig.NumSubframes, regenerate pusch reference
        % resources and carrier with the correct number of slots
        if nSlots > waveConfig.NumSubframes*ofdmInfo.SlotsPerSubframe
            nSubframes = ceil(nSlots/ofdmInfo.SlotsPerSubframe);
            waveConfig.NumSubframes = nSubframes;
            [~,winfo] = nrWaveformGenerator(waveConfig);
            waveformResources = winfo.WaveformResources;
            puschWaveCfg = [];
            puschDefs = [waveConfig.PUSCH{loc}];
            ri = waveformResources.PUSCH(loc);
            ri.PUSCH = puschDefs;
            puschWaveCfg = [puschWaveCfg ri];
            if isempty(puschWaveCfg.PUSCH.NID)
                puschWaveCfg.PUSCH.NID = carrier.NCellID;
            end
        end

        if nSlots == 0
            warning('No scheduled slots found for EVM processing');
            evmInfo(bwpIdx).SubcarrierRMS = [];
            evmInfo(bwpIdx).SubcarrierPeak = [];
            evmInfo(bwpIdx).SymbolRMS = [];
            evmInfo(bwpIdx).SymbolPeak = [];
            evmInfo(bwpIdx).SlotRMS = [];
            evmInfo(bwpIdx).SlotPeak = [];
            evmInfo(bwpIdx).EVMGrid = [];
            evmInfo(bwpIdx).OverallEVM.EV = [];
            evmInfo(bwpIdx).OverallEVM.RMS = [];
            evmInfo(bwpIdx).OverallEVM.Peak = [];
            evmInfo(bwpIdx).Emissions = [];
            evmInfo(bwpIdx).BandwidthPartID = [];
            evmInfo(bwpIdx).AverageFO = 0;
            evmInfo(bwpIdx).AmpImbalance = 0;
            evmInfo(bwpIdx).PhImbalance = 0;
            continue;
        end

        % Array of unallocated PUSCH RBs
        % Each RB is expressed as an offset relative to the
        % nearest edge of the PUSCH allocation
        prbSet = puschWaveCfg.PUSCH.PRBSet;
        if iscolumn(prbSet)
            prbSet = prbSet.';
        end
        nDeltaRB = bwpLen - size(prbSet, 2);
        rbs = (0:carrier.NSizeGrid-1).';
        rbs(prbSet+1)=[];
        rbs(rbs<prbSet(1, 1)) = ...
            rbs(rbs<prbSet(1, 1))- prbSet(1, 1);
        rbs(rbs>prbSet(end, 1)) = ...
            rbs(rbs>prbSet(end, 1))-prbSet(end, 1);
        emissions.DeltaRB = rbs;

        % Generate a reference grid, refGrid, for slots corresponding to
        % the length of the input waveform. This grid contains only the
        % DM-RS and is primarily used for channel estimation
        [refGrid,idealGrid] = hReferenceGrid(carrier,bwpCfg,puschWaveCfg,nSlots,'PUSCH');

        % Slot allocation list for the PUSCH configuration
        activeSlots = [puschWaveCfg.Resources.NSlot];

        % number of FFT Locations in each CP, based on EVM mode (Standard / 3GPP)
        nEVMWindowLocations = 1;
        if evm3GPP
            nEVMWindowLocations = 2;
        end

        % Resize refGrid based on BWP dimensions
        refGrid = refGrid(12*bwpStart+1:12*(bwpStart+bwpLen),:,:);

        % Declare storage variables
        emissions.Absolute = zeros(nDeltaRB, nSlots);   % Array of Absolute emissions , RB vs slot
        emissions.Relative = zeros(nDeltaRB, nSlots);   % Array of Relative emissions , RB vs slot
        rxGridInBandEmissions = [];                     % Demodulated OFDM grid, for In-band emissions
        frameEVM = repmat(hRawEVM([]), 1, max(nFrames,1));   % Per-frame EVM

        % Restrict CP length as per TS 38.101-1 (FR1), TS 38.101-1 (FR2), Annex F.4
        cpLength = double(ofdmInfo.CyclicPrefixLengths(2));

        % Populate puschObj of type 'nrPUSCHConfig'
        % It is to be used for Common Phase error (CPE) estimation
        puschObj = extractPUSCHCfg(waveConfig.PUSCH(loc));
        puschObj = puschObj{1};
        puschObj.NSizeBWP = waveConfig.BandwidthParts{bwpIdx}.NSizeBWP;
        puschObj.NStartBWP = waveConfig.BandwidthParts{bwpIdx}.NStartBWP;

        if isempty(puschObj.NID)
           puschObj.NID =  carrier.NCellID;
        end

        % When evm3GPP is true, two EVM window locations and two CP fractions
        % are selected for 3GPP EVM for OFDM demodulation. If false, a single
        % EVM window location is used, which is centred in the middle of the CP
        if evm3GPP
            W = getEVMWindow(carrier,waveConfig.FrequencyRange,waveConfig.ChannelBandwidth,ofdmInfo.Nfft);
            cpFraction = [0 ; W/cpLength];
        else
            cpFraction = 0.5;      % Use default value
        end

        averageFO = foffsetEstCoarse+integerCfoHz;
        if correctFineFO
            transformPrecoding = puschObj.TransformPrecoding;
            if transformPrecoding
                [rxWaveform,averageFineFO] = hNRFrequencyOffset('fineFO',carrier,bwpCfg,rxWaveform,refGrid,idealGrid,...
                    puschWaveCfg.CDMLengths,sampleRate,activeSlots,puschWaveCfg.Resources,length(puschObj.PRBSet));
            else
                [rxWaveform,averageFineFO] = hNRFrequencyOffset('fineFO',carrier,bwpCfg,rxWaveform,refGrid,idealGrid,...
                    puschWaveCfg.CDMLengths,sampleRate,activeSlots);
            end
            averageFO = averageFO + averageFineFO;
        end

        % OFDM demodulate each slot
        carrier.NSlot = initialNSlot;
        rxGridLow = nrOFDMDemodulate(carrier, rxWaveform, 'CyclicPrefixFraction',cpFraction(1),'SampleRate',sampleRate);
        if nEVMWindowLocations == 2
            rxGridHigh = nrOFDMDemodulate(carrier, rxWaveform, 'CyclicPrefixFraction',cpFraction(2),'SampleRate',sampleRate);
        end
        if (~isempty(emissions.DeltaRB))
            rxGridInBandEmissions = nrOFDMDemodulate(carrier, rxWaveform, 'CyclicPrefixFraction',0.5,'SampleRate',sampleRate);
        end

        % Compute the PUSCH EVM for each active/valid PUSCH slot, store the results
        % in a cell-array for later processing. Skip slots which are not allocated
        slotRange = activeSlots(activeSlots < initialNSlot+nSlots);
        slotRange = slotRange(slotRange >= initialNSlot);
        if isempty(slotRange)
            slotRange = [];
            plotEVM = false;
            warning('No scheduled slots found for EVM processing.');
        end

        % Work only on the BWP on the waveform to simplify indexing
         rxGridLow = rxGridLow(12*bwpStart+1:12*(bwpStart+bwpLen),:,:);
         if evm3GPP
             rxGridHigh = rxGridHigh(12*bwpStart+1:12*(bwpStart+bwpLen),:,:);
         end

        dlFlag = false;
        rxGrids = cell(1,1+evm3GPP);
        rxGrids{1} = rxGridLow;
        if evm3GPP
            rxGrids{2} = rxGridHigh;
        end

        % For each slot, estimate the channel coefficients
        [HestLow,HestHigh]= hChannelEstEVM(rxGrids,refGrid,puschWaveCfg.CDMLengths,nSlots,L,dlFlag);

        bwp = waveConfig.BandwidthParts;
        carrier.NSlot = initialNSlot;

        hest  = cell(1,1+evm3GPP);
        hest{1} = HestLow;
        % create identity here
        % disable_eq = false;
        % if(disable_eq)
        %     hest_id = zeros(size(hest{1}));
        %     [N_SC,N_frames,N_ant,N_port] = size(hest_id);
        %     for i=1:1:N_SC
        %         for j=1:1:N_frames
        %             hest_id(i,j,:,:) = eye(4)+0.0001*1j*eye(4);
        %         end
        %     end
        %     hest{1} = hest_id;
        % end

        if evm3GPP
            rxGrids{2} = rxGridHigh;
            hest{2} = HestHigh;
        end
        [eqSlotGrid,refSlotGrid,emissions,eqSym(:,bwpIdx),refSym(:,bwpIdx)] = hDecodeSlots(carrier,rxGrids,hest,...
            puschWaveCfg,puschObj,bwp{bwpIdx},emissions,rxGridInBandEmissions);

        % Evaluate detailed EVM statistics for a grid of equalized and reference slots
        evmInfoBWP = hEVM(carrier,eqSlotGrid,refSlotGrid);
        evm = evmInfoBWP.EVM;

        % Add emissions to evmInfo
        evmInfoBWP.Emissions = emissions;

        evmInfoBWP.BandwidthPartID = waveConfig.BandwidthParts{bwpIdx}.BandwidthPartID;
        % Remove EVM field as it not part of the output
        evmInfoBWP = rmfield(evmInfoBWP,'EVM');
        evmInfoBWP.AverageFO = averageFO;
        evmInfoBWP.AmpImbalance = ampImbEst;
        evmInfoBWP.PhImbalance = phImbEst;
        evmInfo = [evmInfo;evmInfoBWP];

        % Display per slot, and per EVM edge EVM statistics
        if displayEVM
            disp("EVM stats for BWP idx : " + num2str(bwpIdx));
        end
        for slotIdx  = slotRange
            for e = 1:nEVMWindowLocations
                if (e == 1)
                    edge = 'Low edge';
                    if nEVMWindowLocations == 1
                        edge = '';              % Print only single EVM per slot
                    end
                else
                    edge = 'High edge';
                end
                if displayEVM
                    fprintf('%s RMS EVM, Peak EVM, slot %d: %0.3f %0.3f%%\n',edge,slotIdx,evm(e, slotIdx+1).RMS*100,evm(e, slotIdx+1).Peak*100);
                end
            end
        end

        printFrameAvg = 1;      % Ensures only fully occupied frames are printed
        % After we've filled a frame or if we're at the end of a signal
        % shorter than a frame, do EVM averaging
        if (nFrames == 0)
            nFrames = 1;     % Below loop needs to run at least once
            printFrameAvg = 0; % Don't print low/high EVM as we don't have sufficient slots to fill up a frame
        end

        % 1-based indexing for accessing evm
        % Limit frame-averaging to complete frames only
        slotRange = slotRange+1;
        slotRange = slotRange(slotRange <= (nFrames*waveConfig.NumSubframes*ofdmInfo.SlotsPerSubframe));

        % loop through each frame, selecting the frames with higher RMS( when
        % measuring 3GPP EVM)
        for frameIdx = 0:nFrames-1
            frameLowEVM = hRawEVM(cat(1,evm(1,slotRange).EV));
            frameEVM(frameIdx+1) = frameLowEVM;
            if evm3GPP
                frameHighEVM = hRawEVM(cat(1,evm(2,slotRange).EV));
                if frameHighEVM.RMS > frameLowEVM.RMS
                    frameEVM(frameIdx+1) = frameHighEVM;
                end
            end
            if printFrameAvg && displayEVM
                if evm3GPP
                    fprintf('Averaged low edge RMS EVM, frame %d: %0.3f%%\n',frameIdx,frameLowEVM.RMS*100);
                    fprintf('Averaged high edge RMS EVM, frame %d: %0.3f%%\n',frameIdx,frameHighEVM.RMS*100);
                end
                fprintf('Averaged RMS EVM frame %d: %0.3f%%\n',frameIdx,frameEVM(frameIdx+1).RMS*100);
            end
        end
        if displayEVM
            fprintf('Averaged overall RMS EVM: %0.3f%%\n', evmInfo(bwpIdx).OverallEVM.RMS*100);
        end
        if plotEVM
            % Update the plots (symbol,subcarrier,slot and grid-wise)
            hEVMPlots(evmInfo(bwpIdx),eqSym{1,bwpIdx},refSym{1,bwpIdx},'PUSCH');
        end
        if displayEVM
            disp("Peak EVM : " + string((evmInfo(bwpIdx).OverallEVM.Peak)*100) + "%");
        end
    end
end

function W = getEVMWindow(carrier,frequencyRange,channelBandwidth,nFFT)
%   W = getEVMWindow(CARRIER,FREQUENCYRANGE,CHANNELBANDWIDTH,NFFT) is the
%   error vector magnitude window length, as mentioned in TS 38.101-1(2), Section
%   F.5 . W is defined for a given combination of subcarrier
%   spacing, channel bandwidth/fft length, frequency range and CP type. W
%   is subsequently used as an intermediate value to decide the CP Fraction
%   for OFDM demodulation

    scsFR1 = [15 30 60];
    scsFR2 = [60 120];
    % BW MHz        5  10 15 20  25  30  40  50  60  70  80  90  100
    nfftFR1      = [256 384 512 768 1024 1536 2048 3072  4096];
    WsFR1 =          [NaN   NaN    18   NaN     36     54     72    108     144; 	%NormalCp,15kHz
                        9   NaN    18    27     36     54     72    108     144;	%NormalCp,30kHz
                        9    14    18    27     36     54     72    NaN     NaN;    %NormalCp,60kHz
                       54    80   106   164    220    330    440    NaN     NaN];   %ExtendedCp,60kHz

    % BW MHz        50 100 200 400
    nfftFR2      = [512 1024 2048 4096];
    WsFR2 = [  NaN  36  72  144;		 %NormalCP,60kHz
                18  36  72  144;    	 %NormalCP,120kHz
               NaN 220 440  880];        %ExtendedCP,60kHz
    W = []; %#ok<NASGU>
    if (strcmpi(frequencyRange,'FR1'))
        rowIdx = find(carrier.SubcarrierSpacing == scsFR1) + double(strcmpi(carrier.CyclicPrefix,'extended'));
        W = WsFR1(rowIdx,nFFT == nfftFR1);
        if channelBandwidth == 25
            if nFFT == 512 && carrier.SubcarrierSpacing == 60
                if strcmpi(carrier.CyclicPrefix,'extended')
                    W = 110;
                else
                    W = 18;
                end
            end
        end
    else
        rowIdx = find(carrier.SubcarrierSpacing == scsFR2) + double(strcmpi(carrier.CyclicPrefix,'extended'))*2;
        W = WsFR2(rowIdx,nFFT == nfftFR2);
    end
    % Filter out invalid combinations
    if isnan(W) || isempty(W)
        error('Invalid FFT/SCS/BW combination');
    end
end

function [puschObj] = extractPUSCHCfg(puschWaveCfg)
    % Extract relevant parameters from waveform generator PUSCH struct to build an object of type 'nrPUSCHConfig' 
    puschObj = cell(1);
    for idx = 1:length(puschWaveCfg)
        puschObj{idx} = nrPUSCHConfig;
        puschObj{idx}.Modulation                     = puschWaveCfg{idx}.Modulation;
        puschObj{idx}.NumLayers                      = puschWaveCfg{idx}.NumLayers;
        puschObj{idx}.MappingType                    = puschWaveCfg{idx}.MappingType;
        puschObj{idx}.SymbolAllocation               = puschWaveCfg{idx}.SymbolAllocation;
        puschObj{idx}.PRBSet                         = puschWaveCfg{idx}.PRBSet;
        puschObj{idx}.TransformPrecoding             = puschWaveCfg{idx}.TransformPrecoding;
        puschObj{idx}.TransmissionScheme             = puschWaveCfg{idx}.TransmissionScheme;
        puschObj{idx}.NumAntennaPorts                = puschWaveCfg{idx}.NumAntennaPorts;
        puschObj{idx}.TPMI                           = puschWaveCfg{idx}.TPMI;
        puschObj{idx}.FrequencyHopping               = puschWaveCfg{idx}.FrequencyHopping;
        puschObj{idx}.SecondHopStartPRB              = puschWaveCfg{idx}.SecondHopStartPRB;
        puschObj{idx}.NID                            = puschWaveCfg{idx}.NID;
        puschObj{idx}.RNTI                           = puschWaveCfg{idx}.RNTI;

        % Set DM-RS parameters
        puschObj{idx}.DMRS                           = puschWaveCfg{idx}.DMRS;

        % Set PT-RS parameters
        puschObj{idx}.EnablePTRS                     = puschWaveCfg{idx}.EnablePTRS;
        puschObj{idx}.PTRS                           = puschWaveCfg{idx}.PTRS;
    end
end

function cfgObj = mapStruct2ObjUplink(cfgS)

    % Map the configurations into object-based ones needed by nrWaveformGenerator
    cfgObj = nrULCarrierConfig;

    % Top-level parameters
    [cfgObj,cfgS] = mapStruct2ObjCCCommon(cfgObj,cfgS);

    % Continue with uplink channels and resources 
    cfgObj.PUSCH                 = mapPUSCH(cfgS);

end

% Component carrier common parameters
function [cfgObj,cfgS] = mapStruct2ObjCCCommon(cfgObj,cfgS)

    % Top-level parameters
    cfgObj.FrequencyRange        = cfgS.FrequencyRange;
    cfgObj.ChannelBandwidth      = cfgS.ChannelBandwidth;
    cfgObj.NCellID               = cfgS.NCellID;
    cfgObj.NumSubframes          = cfgS.NumSubframes;
    if isprop(cfgObj,'WindowingPercent') && isfield(cfgS,'Windowing')
        % Only map to either off or default windowing levels
        if cfgS.Windowing == 0
            w = 0;
        else
            w = [];
        end
        cfgObj.WindowingPercent = w;
    end
    % Add in a name field should one not exist and copy across string otherwise
    if ~isfield(cfgS,'Name')
        cfgS.Name = '';
    else
        cfgObj.Label = cfgS.Name;
    end

    cfgObj.SCSCarriers           = mapCarriers(cfgS);
    [cfgObj.BandwidthParts,cfgS] = mapBWP(cfgS);

end

function scsCfg = mapCarriers(cfgS)

    % SCS carrier configuration
    carriers = cfgS.Carriers;
    for idx = 1:length(carriers)
        scsCfg{idx} = nrSCSCarrierConfig;
        % Map across parameters
        scsCfg{idx}.SubcarrierSpacing = carriers(idx).SubcarrierSpacing; %#ok<*AGROW>
        scsCfg{idx}.NSizeGrid         = carriers(idx).NRB;
        scsCfg{idx}.NStartGrid        = carriers(idx).RBStart;
    end
end

function [bwp2,cfgS] = mapBWP(cfgS)
  
    bwp2 = {};

    bwp = cfgS.BWP;
    for idx = 1:length(bwp)
        bwp2{idx} = nrWavegenBWPConfig;
        bwp2{idx}.BandwidthPartID = idx;
        bwp2{idx}.SubcarrierSpacing = bwp(idx).SubcarrierSpacing;
        bwp2{idx}.CyclicPrefix      = bwp(idx).CyclicPrefix;
        bwp2{idx}.NSizeBWP          = bwp(idx).NRB;
        % Find carrier for corresponding bwp 
        carrierid = [cfgS.Carriers(:).SubcarrierSpacing] == bwp2{idx}.SubcarrierSpacing;
        bwp2{idx}.NStartBWP         = bwp(idx).RBOffset + cfgS.Carriers(carrierid).RBStart;
        cfgS.BWP(idx).Carrier = find(carrierid,1);     % Link back the SCS carrier index associated with the BWP (linked through SCS value)
    end
end

function pusch2 = mapPUSCH(cfgS)
    pusch2 = {};
    % PUSCH
    pusch = cfgS.PUSCH;
    for idx = 1:length(pusch)
        pusch2{idx} = nrWavegenPUSCHConfig;
        % General shared parameters
        pusch2{idx}.Enable               = pusch(idx).Enable;
        pusch2{idx}.Label                = pusch(idx).Name;
        pusch2{idx}.Power                = pusch(idx).Power;
        pusch2{idx}.BandwidthPartID      = pusch(idx).BWP;
        pusch2{idx}.Coding               = pusch(idx).EnableCoding;
        pusch2{idx}.DataSource           = pusch(idx).DataSource;
        pusch2{idx}.TargetCodeRate       = pusch(idx).TargetCodeRate;
        pusch2{idx}.XOverhead            = pusch(idx).Xoh_PUSCH;
        pusch2{idx}.Modulation           = pusch(idx).Modulation;
        pusch2{idx}.NumLayers            = pusch(idx).NLayers;
        pusch2{idx}.RVSequence           = pusch(idx).RVSequence;
        pusch2{idx}.SymbolAllocation     = [pusch(idx).AllocatedSymbols(1) (pusch(idx).AllocatedSymbols(end)+1-pusch(idx).AllocatedSymbols(1))];
        pusch2{idx}.SlotAllocation       = pusch(idx).AllocatedSlots;
        pusch2{idx}.Period               = pusch(idx).AllocatedPeriod;
        pusch2{idx}.PRBSet               = pusch(idx).AllocatedPRB;
        pusch2{idx}.RNTI                 = pusch(idx).RNTI;
        pusch2{idx}.NID                  = pusch(idx).NID;
        pusch2{idx}.MappingType          = pusch(idx).PUSCHMappingType;
  
        % PUSCH specific parameters
        pusch2{idx}.TransformPrecoding   = pusch(idx).TransformPrecoding;
        pusch2{idx}.TransmissionScheme   = pusch(idx).TxScheme;
        pusch2{idx}.NumAntennaPorts      = pusch(idx).NAntennaPorts;
        pusch2{idx}.TPMI                 = pusch(idx).TPMI;
        if strcmpi(pusch(idx).InterSlotFreqHopping,'enabled')
            pusch2{idx}.FrequencyHopping = 'interSlot';
        elseif strcmpi(pusch(idx).IntraSlotFreqHopping,'enabled')
            pusch2{idx}.FrequencyHopping = 'intraSlot';
        else
            pusch2{idx}.FrequencyHopping = 'neither';
        end
        if ~strcmpi(pusch2{idx}.FrequencyHopping,'neither')
            if isempty(pusch(idx).RBOffset) % Replace the empty value with 0
                pusch(idx).RBOffset = 0;
            end
            secondHopPRB = mod(pusch(idx).AllocatedPRB(1)+pusch(idx).RBOffset,cfgS.BWP(idx).NRB);
            pusch(idx).RBOffset = secondHopPRB; % Replace RB offset with the calculated secondHopPRB
        end
        pusch2{idx}.SecondHopStartPRB    = pusch(idx).RBOffset;

        % UCI part
        pusch2{idx}.EnableULSCH          = ~pusch(idx).DisableULSCH;

        % Antenna port and DM-RS configuration (TS 38.211 section 7.4.1.1)
        % Using DM-RS reference point CRB0
        dmrs = nrPUSCHDMRSConfig;
        dmrs.NRSID                     = pusch(idx).NRSID;
        dmrs.NumCDMGroupsWithoutData   = pusch(idx).NumCDMGroupsWithoutData;
        dmrs.DMRSTypeAPosition         = pusch(idx).DMRSTypeAPosition;
        dmrs.DMRSAdditionalPosition    = pusch(idx).DMRSAdditionalPosition;
        dmrs.DMRSLength                = pusch(idx).DMRSLength;
        dmrs.DMRSPortSet               = pusch(idx).PortSet;
        dmrs.NIDNSCID                  = pusch(idx).NIDNSCID;
        dmrs.NSCID                     = pusch(idx).NSCID;
        dmrs.DMRSConfigurationType     = pusch(idx).DMRSConfigurationType;

        if isfield(pusch(idx),'GroupHopping') && pusch(idx).TransformPrecoding
            if strcmpi(pusch(idx).GroupHopping,'enable')
                dmrs.GroupHopping = 1;
            elseif strcmpi(pusch(idx).GroupHopping,'disable')
                dmrs.SequenceHopping = 1;
            end
        end
        if isfield(pusch(idx),'NumCDMGroupsWithoutData') && ~isempty(pusch(idx).NumCDMGroupsWithoutData) ...
                && ~pusch(idx).TransformPrecoding && pusch(idx).NumCDMGroupsWithoutData
            dmrs.NumCDMGroupsWithoutData = pusch(idx).NumCDMGroupsWithoutData;
        else
            dmrs.NumCDMGroupsWithoutData = 1 + pusch(idx).TransformPrecoding;
        end
        pusch2{idx}.DMRSPower             = pusch(idx).PowerDMRS;
        pusch2{idx}.DMRS = dmrs;

        numDMRSPorts = numel(dmrs.DMRSPortSet);
        if numDMRSPorts
            % Assign the number of DM-RS antenna ports to layers
            pusch2{idx}.NumLayers = numDMRSPorts;
        end

        % PT-RS configuration
        pusch2{idx}.EnablePTRS              = pusch(idx).EnablePTRS;
        ptrs{idx} = nrPUSCHPTRSConfig;
        ptrs{idx}.TimeDensity               = pusch(idx).PTRSTimeDensity;
        ptrs{idx}.FrequencyDensity          = pusch(idx).PTRSFrequencyDensity;
        ptrs{idx}.REOffset                  = pusch(idx).PTRSREOffset;
        ptrs{idx}.PTRSPortSet               = pusch(idx).PTRSPortSet;
        ptrs{idx}.NID                       = pusch(idx).PTRSNID;
        pusch2{idx}.PTRSPower               = pusch(idx).PowerPTRS;
        pusch2{idx}.PTRS = ptrs{idx};
    end
end

function validateInputs(evmCfg)
%   validateInputs(EVMCFG)
%   validates the evm configuration
%   EVMCFG used for EVM measurement.

    fcnName = 'hNRPUSCHEVM';

    % Validate 'evmCfg'
    if isfield(evmCfg,'Evm3GPP')
        validateattributes(evmCfg.Evm3GPP,{'logical','double'},{'nonempty'},fcnName,'Evm3GPP');
    end
    if isfield(evmCfg,'PlotEVM')
        validateattributes(evmCfg.PlotEVM,{'logical','double'},{'nonempty'},fcnName,'PlotEVM');
    end
    if isfield(evmCfg,'DisplayEVM')
        validateattributes(evmCfg.DisplayEVM,{'logical','double'},{'nonempty'},fcnName,'DisplayEVM');
    end
    if isfield(evmCfg,'InitialNSlot')
        validateattributes(evmCfg.InitialNSlot,{'numeric'},{'nonnegative'},fcnName,'InitialNSlot');
    end
    if isfield(evmCfg,'SampleRate')
        validateattributes(evmCfg.SampleRate,{'numeric'},{'integer','positive'},fcnName,'SampleRate');
    end
    if isfield(evmCfg,'CorrectCoarseFO')
        validateattributes(evmCfg.CorrectCoarseFO,{'logical','double'},{'nonempty'},fcnName,'CorrectCoarseFO');
    end
    if isfield(evmCfg,'CorrectFineFO')
        validateattributes(evmCfg.CorrectFineFO,{'logical','double'},{'nonempty'},fcnName,'CorrectFineFO');
    end
    if isfield(evmCfg,'TimeSyncEnable')
        validateattributes(evmCfg.TimeSyncEnable,{'logical','double'},{'nonempty'},fcnName,'TimeSyncEnable');
    end
    if isfield(evmCfg,'UseWholeGrid')
        validateattributes(evmCfg.UseWholeGrid,{'logical','double'},{'nonempty'},fcnName,'UseWholeGrid');
    end
    if isfield(evmCfg,'IQImbalance')
        validateattributes(evmCfg.IQImbalance,{'logical','double'},{'nonempty'},fcnName,'IQImbalance');
    end
end