clc
clear all
L=10
snr_dB = -10:1:10;
N = 64;
L = 10;
Tu = 64;
Tg = Tu/4;

BPSK = [1 -1];
BER_bpsk = zeros(1,length(snr_dB));
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


    %Measure number of bits in each frame
    num_bits = length(bits);


    %Perform BPSK Modulation
    bpsk_mod = [];
    bpsk_mod = 2*bits -1;


    %Convert Serial data to parallel data
    parallel_data = reshape(bpsk_mod, [], 64);

    %Convert data to time domain
    ofdmSymbols_bpsk = ifft(parallel_data,[],2);
    
    %Add Cyclic Prefix
    size_bpsk = size(ofdmSymbols_bpsk);
    ofdm_bpsk_cp = zeros(size_bpsk(1),74);
    for num = 1:size(ofdmSymbols_bpsk, 1)
        data = ofdmSymbols_bpsk(num, :);
        cyclic_prefix = data (end-10+1:end);
        data_with_prefix= [cyclic_prefix data];
        ofdm_bpsk_cp(num, :) = data_with_prefix;
    end
    
    %Create the Channel
    
    BPSK = unique(bpsk_mod);
    BPSK_SNR_db = snr_dB(loop);
    BPSK_symbols = BPSK;
    % calculate the magnitude of each complex number
    BPSK_magnitudeArray = abs(BPSK_symbols);
    % square the magnitude of each complex number and sum them
    BPSK_sumOfSquares = sum(BPSK_magnitudeArray.^2);
    % divide the sum by the length of the original array
    BPSK_c = BPSK_sumOfSquares / length(BPSK_symbols);
    BPSK_SNR_W = db2pow(BPSK_SNR_db);
    BPSK_var_w = BPSK_c/BPSK_SNR_W;
    bpsk_txSig = ofdm_bpsk_cp;
    taps=10;  %number of channel taps
    bpsk_channel = [randn(1,taps)+j*randn(1,taps) zeros(1,size(bpsk_txSig,2)-taps)].' ; % The Rayleigh fading channel padded with zeros
    bpsk_channel_mat = toeplitz(bpsk_channel); %  Convolution is the same as multiplying with the toeplitz matrix.
    
    %Pass the OFDM symbols through the channel
    bpsk_rcvsig=bpsk_channel_mat*bpsk_txSig'; % Pass through the Rayleigh fading channel
    bpsk_sig = bpsk_rcvsig + sqrt(BPSK_var_w/2)*(randn(size(bpsk_rcvsig)) + 1i*randn(size(bpsk_rcvsig)));
    bpsk_sig_equal = inv(bpsk_channel_mat'*bpsk_channel_mat)*(bpsk_channel_mat')*bpsk_sig; % 0  channel normalization
    y_bpsk = bpsk_sig_equal';
     
    %Remove Cyclic prefix from received signal
    bpsk_received_cp_removed = y_bpsk(:,L+ 1:end);

    %Convert signal to frequency domain
    bpsk_received = fft(bpsk_received_cp_removed, [], 2);

    %Convert signal from parallel to serial data
    serial_data = reshape(bpsk_received, [], 1);

    %Demodulate the signal
   BPSK_bits_Rx = [];

   % Determine if positive or negative
   for index = 1:length(serial_data)
       if serial_data(index) > 0
           BPSK_bits_Rx = [BPSK_bits_Rx 1];
       elseif serial_data(index) < 0
           BPSK_bits_Rx = [BPSK_bits_Rx 0];
       end
   end

    %Recreate video
         orig_class = class(currentFrame);
     orig_size = size(currentFrame);
     bpsk_reconstructedFrames(:,:,:,i) = reshape(typecast(uint8(bin2dec ...
         (char(reshape(BPSK_bits_Rx, 8, [])+'0').')),orig_class), orig_size);

     % Calculate BER
     BER_bpsk(loop) = BER_bpsk(loop) + sum(BPSK_bits_Rx ~= bits.');
    %disp(i)    
end

% Normalize the BER for every frame. 
BER_bpsk(loop) = BER_bpsk(loop)/(num_bits*size_frame);
disp(loop)

end 

%plot BER vs. SNR
figure(1)
semilogy(snr_dB, BER_bpsk);
xlabel('SNR (dB)')
ylabel('BPSK BER (Pe)')
title('BPSK Probability of error vs. SNR')
savefig('BPSK_BER.fig')

% play the video
implay(bpsk_reconstructedFrames);