% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren

clear

% Initialization
EbN0_db = 0:10;                     % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 20;             % Size of training sequence (in nr bits)
nr_blocks = 100;                     % The number of blocks to simulate
Q = 8;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);
%pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(1, length(EbN0_db));   % Error counter

nr_sync_errors = zeros(1, length(EbN0_db)); %syncronization error
square_phase_error = zeros(1, length(EbN0_db)); %phase errors
forced_errors=0;
nr_errors_forced = zeros(1, length(forced_errors));   % Error counter

for snr_point = 1:length(EbN0_db)  
  for k = 1:length(nr_training_bits)
      % Loop over several blocks to get sufficient statistics.
      for blk = 1:nr_blocks

        %%%
        %%% Transmitter
        %%%

        % Generate training sequence.
        b_train = training_sequence(nr_training_bits(k));        
        
        % Generate random source data {0, 1}.
        b_data = random_data(nr_data_bits);

        % Generate guard sequence.
        b_guard = random_data(nr_guard_bits);

        % Multiplex training and data into one sequence.
        b = [b_guard b_train b_data b_guard];
        % Map bits into complex-valued QPSK symbols.
        d = qpsk(b);

        % Upsample the signal, apply pulse shaping.
        tx = upfirdn(d, pulse_shape, Q, 1);

        %%%
        %%% AWGN Channel
        %%%

        % Compute variance of complex noise according to report.
        sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);

        % Create noise vector.
        n = sqrt(sigma_sqr/2)*(randn(size(tx))+1i*randn(size(tx)));

        % Received signal.
        rx = tx + n;

        %%%
        %%% Receiver
        %%%

        % Matched filtering.
        mf=conv(mf_pulse_shape,rx);

        % Down sampling. t_samp is the first sample, the remaining samples are all

        % Synchronization. The position and size of the search window
        % is here set arbitrarily. Note that you might need to change these
        % parameters. Use sensible values (hint: plot the correlation
        % function used for syncing)! 
        t_start=Q*(nr_guard_bits-2)/2;
        t_end=t_start+50;
        t_samp = sync(mf, b_train, Q, t_start, t_end);
        %t_samp = 1+Q*nr_guard_bits/2+Q-1; %Force the value of t_samp, we can add error here
        nr_sync_errors(snr_point) = nr_sync_errors(snr_point) + abs(t_samp-(1+Q*nr_guard_bits/2+Q-1));

        % separated by a factor of Q. Only training+data samples are kept.
        r = mf(t_samp:Q:t_samp+Q*(nr_training_bits(k)+nr_data_bits)/2-1);


        % Phase estimation and correction.
        phihat = phase_estimation(r, b_train);
        %phihat=0; %Force the value of phihat, we can add error here
        square_phase_error(snr_point) = square_phase_error(snr_point)+abs(phihat)^2; 
        r = r * exp(-1i*phihat);
        
        %Plot constellation
        %figure()
        %plot(r);

        % Make decisions. Note that dhat will include training sequence bits
        % as well.
        bhat = detect(r);

        % Count errors. Note that only the data bits and not the training bits
        % are included in the comparison. The last data bits are missing as well
        % since the whole impulse response due to the last symbol is not
        % included in the simulation program above.
        temp=bhat(1+nr_training_bits(k):nr_training_bits(k)+nr_data_bits) ~= b_data;
        nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);
        %nr_errors_forced(k) = nr_errors_forced(k) + sum(temp);

        % Next block.
      end
  end
  % Next Eb/No value.
end

% Compute the BER. 
BER = nr_errors / nr_data_bits / nr_blocks;
SyncErrorRate = nr_sync_errors/ nr_blocks;
square_phase_error = sqrt(square_phase_error/nr_blocks);
BER_forced = nr_errors_forced / nr_data_bits / nr_blocks;

%Display BER infos
% disp(BER);
% figure(2)
% plot(EbN0_db,BER,'b')
% hold on;
% plot(EbN0_db,1-normcdf(sqrt(2*10.^(EbN0_db/10))),'r');
% hold off;

%Display Sync error infos
% disp(SyncErrorRate);
% plot(EbN0_db,SyncErrorRate,'b');


%Display phase error infos
%disp(square_phase_error)
%plot(EbN0_db,square_phase_error);

%Display BER under added errors
%disp(BER_forced)
%plot(forced_errors,BER_forced);


%Display influence of training bits
%plot(nr_training_bits,SyncErrorRate);
%plot(nr_training_bits,square_phase_error);

