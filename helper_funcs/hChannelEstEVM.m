function [HestLow,HestHigh] = hChannelEstEVM(rxGrids,refGrid,cdmLengths,nSlots10Ms,L,dlFlag)
    %  [HESTLOW,HESTHIGH] = hChannelEstEVM(RXGRIDS,REFGRID,CDMLENGTHS,NSLOTS10MS,L,DLFLAG)
    %  Estimates the channel returning the channel coefficients HESTLOW.
    %  HESTLOW is a K-by-N-by-R-by-P array where K is the number of
    %  subcarriers, N is the number of symbols, and R is the number of
    %  receive antennas and P is the number of reference signal ports.
    %  HESTHIGH has same dimensions as HestLow and is additionally returned
    %  when EVM3GPP is true.
    %  RXGRIDS is a cell array of demodulated IQs of length 1 or 2
    %  dependant on the mode of channel estimation, 3GPP or non-3GPP. Each
    %  element of RXGRIDS must be an array of size K-by-N-by-R.
    %  REFGRID is a predefined reference array with nonzero elements
    %  representing the reference symbols in their appropriate locations.
    %  It is of size K-by-N-by-P. REFGRID can span multiple slots.
    %  CDMLENGTHS is a 2-element row vector [FD TD] specifying the length
    %  of FD-CDM and TD-CDM despreading to perform.
    %  NSLOTS10MS is the number of slots in a span of 10ms
    %  L is the number of symbols in a slot.
    %  DLFLAG when set to true, enables smoothening of the channel
    %  coefficients in the frequency direction, using a moving average
    %  filter. Smoothening is performed as described in TS 38.104 Annex B.6
    %  (FR1) or C.6 (FR2). When this parameter is not specified, the time
    %  averaging across the duration of REFGRID is enabled.
    
    % Copyright 2021 The MathWorks, Inc.
    
    HestLow = [];                       % Channel estimation for 1st CP position
    HestHigh = [];                      % Channel estimation for 2nd CP position
    if isempty(refGrid)
        return;
    end

    % Extract channel estimation mode
    evm3GPP = length(rxGrids)>1;

    % Extract rxGridLow and if required, rxGridHigh
    rxGridLow = rxGrids{1};
    if evm3GPP
        rxGridHigh = rxGrids{2};
    end

    nSlots = floor(size(refGrid,2)/L);
    if evm3GPP
        nBlocks = ceil(nSlots/nSlots10Ms);
        nIterations = nSlots;
        
        % For downlink 3GPP channel estimation, the averaging needs to be
        % done in chunks of 10ms. For uplink channel estimation, averaging
        % is done on a per-lsot basis.
        if dlFlag
            nIterations = nBlocks;
        end
    
        % For each 10 ms block, estimate the channel coefficients
        for idx = 1:nIterations
            % Index of symbols within current block. If a symbol index
            % exceeds the length of the reference grid, remove it
            if dlFlag
                symIdx = (idx-1)*L*nSlots10Ms+(1:(L*nSlots10Ms));
            else
                symIdx = (idx-1)*L+1:idx*L;
            end
            symIdx(symIdx>size(refGrid, 2)) = [];

            % Use a smoothing filter in the frequency direction when dlFlag
            % is true
            if dlFlag
                HestLowBlk = hChannelEstimateEVM3GPP(rxGridLow(:, symIdx, :),refGrid(:, symIdx,:),'movingAvgFilter','CDMLengths',cdmLengths);
                HestHighBlk = hChannelEstimateEVM3GPP(rxGridHigh(:, symIdx, :),refGrid(:, symIdx,:),'movingAvgFilter','CDMLengths',cdmLengths);
            else
                HestLowBlk = hChannelEstimateEVM3GPP(rxGridLow(:, symIdx, :),refGrid(:, symIdx,:),'CDMLengths',cdmLengths);
                HestHighBlk = hChannelEstimateEVM3GPP(rxGridHigh(:, symIdx, :),refGrid(:, symIdx,:),'CDMLengths',cdmLengths);
            end
            HestLow  = [HestLow HestLowBlk]; %#ok<*AGROW> 
            HestHigh = [HestHigh HestHighBlk];
        end
    else
        % Compute channel estimates for each slot
        for slotIdx = 1:nSlots
            symIdx = (slotIdx-1)*L+1:slotIdx*L;
            symIdx(symIdx>size(rxGridLow, 2)) = [];
            HestLowSlot = nrChannelEstimate(rxGridLow(:, symIdx, :),refGrid(:, symIdx,:),'CDMLengths',cdmLengths);
            HestLow = [HestLow HestLowSlot];
        end
    end
end