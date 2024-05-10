
clc
clear all


N = 64;
L = 10;
Tu = 64;
Tg = Tu/4;
snr_dB = -10:1:10;
QPSK = [0.7071+0.7071i -0.7071+0.7071i -0.7071-0.7071i 0.7071-0.7071i];
BER_qpsk = zeros(1,length(snr_dB));

% Reading the video 
v = VideoReader('test_video.mp4');
frames = read(v,[1 100]); % reading the first 200 frames of the video
size_frame = size(frames,4);

for loop = 1:1:length(snr_dB)

for i = 1:1:size(frames,4)
    % Considering the current frame
    currentFrame_o = squeeze(frames(:,:,:,i));
    % Reducing the resolution
    currentFrame = imresize(currentFrame_o, [size(frames,1)/4, size(frames,2)/4]);
    % Converting the video frames to bits
    bits = reshape((dec2bin(typecast(squeeze ...
        (currentFrame(:)), 'uint8'), 8) - '0').', 1,[]);

    % Measure number of bits in a frame
    num_bits = length(bits);

    % Perform QPSK Modulation
    qpsk_mod = [];
    bits_per_symb_qpsk = 2;
    
    %Create pairs of bits 
    grouped_qpsk = reshape(bits,bits_per_symb_qpsk,[]).';
  
   % Map bit pairs to QPSK symbols (every pi/4 phase)
    qpsk_mod = exp(1j*(pi/4)*(2*grouped_qpsk(:,1) + grouped_qpsk(:,2)*2-3));

    %Convert serial data to parallel data 
    parallel_data = reshape(qpsk_mod, [], 64);

    % Convert symbols from frequency to time domain
    ofdmSymbols_qpsk = ifft(parallel_data,[],2);
    
    % Add cyclic prefix
    size_qpsk = size(ofdmSymbols_qpsk);
    ofdm_qpsk_cp = zeros(size_qpsk(1),74);
    for num = 1:size(ofdmSymbols_qpsk, 1)
        data = ofdmSymbols_qpsk(num, :);
        cyclic_prefix = data (end-10+1:end);
        data_with_prefix= [cyclic_prefix data];
        ofdm_qpsk_cp(num, :) = data_with_prefix;
    end
    
    
    % Create the channel

    QPSK = unique(qpsk_mod);
    QPSK_SNR_db = snr_dB(loop);
    QPSK_symbols = QPSK;
    % calculate the magnitude of each complex number
    QPSK_magnitudeArray = abs(QPSK_symbols);
    % square the magnitude of each complex number and sum them
    QPSK_sumOfSquares = sum(QPSK_magnitudeArray.^2);
    % divide the sum by the length of the original array
    QPSK_c = QPSK_sumOfSquares / length(QPSK_symbols);
    QPSK_SNR_W = db2pow(QPSK_SNR_db);
    QPSK_var_w = QPSK_c/QPSK_SNR_W;
    qpsk_txSig = ofdm_qpsk_cp;
    taps=10;  %number of channel taps
    qpsk_channel = [randn(1,taps)+j*randn(1,taps) zeros(1,size(qpsk_txSig,2)-taps)].' ; % The Rayleigh fading channel padded with zeros
    qpsk_channel_mat = toeplitz(qpsk_channel); %  Convolution is the same as multiplying with the toeplitz matrix.
    

    %Pass the signal through the channel
    qpsk_rcv_sig=qpsk_channel_mat*qpsk_txSig'; % Pass through the Rayleigh fading channel
    qpsk_sig = qpsk_rcv_sig + sqrt(QPSK_var_w/2)*(randn(size(qpsk_rcv_sig)) + 1i*randn(size(qpsk_rcv_sig)));
    qpsk_sig_equal = inv(qpsk_channel_mat'*qpsk_channel_mat)*(qpsk_channel_mat')*qpsk_sig; % 0 forcing channel normalization
    y_qpsk = qpsk_sig_equal';
     
    %Remove Cyclic Prefix of the received signal
    qpsk_received_cp_removed = y_qpsk(:,L+ 1:end);

    %Convert time domain signal to frequency domain
    qpsk_received = fft(qpsk_received_cp_removed, [], 2);

    %Convert parallel data to serial data
    serial_data = reshape(qpsk_received, [], 1);

    %Demodulate the signal

    QPSK_bits_Rx = [];
   % Determine quadrant
   for index = 1:length(serial_data)
       if real(serial_data(index)) < 0 && imag(serial_data(index)) < 0
           QPSK_bits_Rx = [QPSK_bits_Rx 0 0];
       elseif real(serial_data(index)) > 0 && imag(serial_data(index)) < 0
           QPSK_bits_Rx = [QPSK_bits_Rx 1 0];
       elseif real(serial_data(index)) > 0 && imag(serial_data(index)) > 0
           QPSK_bits_Rx = [QPSK_bits_Rx 1 1];
       elseif real(serial_data(index)) < 0 && imag(serial_data(index)) > 0
           QPSK_bits_Rx = [QPSK_bits_Rx 0 1];
       end
   end

    %Recreate video
         orig_class = class(currentFrame);
     orig_size = size(currentFrame);
     qpsk_reconstructedFrames(:,:,:,i) = reshape(typecast(uint8(bin2dec ...
         (char(reshape(QPSK_bits_Rx, 8, [])+'0').')),orig_class), orig_size);

     %Calculate BER
    BER_qpsk(loop) = BER_qpsk(loop) + sum(QPSK_bits_Rx ~= bits.');
    %disp(i)

end

%Normalize BER for every frame
BER_qpsk(loop) = BER_qpsk(loop)/(num_bits*size_frame);
disp(loop)
end 

%plot BER vs. SNR 

figure(1)
semilogy(snr_dB, BER_qpsk);
xlabel('SNR (dB)');
ylabel('QPSK BER (Pe)');
title('QPSK Probability of error vs. SNR');
savefig('QPSK_BER.fig');


% play the video
implay(qpsk_reconstructedFrames);