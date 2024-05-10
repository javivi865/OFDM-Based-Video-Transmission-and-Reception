clc
clear all
L=10;
N = 64;
Tu = 64;
snr_dB = 5;

QAM = [1+1i -1+1i -1-1i 1-1i 3+3i -3+3i -3-3i 3-3i 3+1i -3+1i -3-1i 3-1i 1+3i -1+3i -1-3i 1-3i];
gray_map = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10]; % Gray mapping table for 16 QAM


BER_qam = zeros(1,length(snr_dB));
BER_qam_3 = zeros(1,length(snr_dB), 'double');
BER_qam_5 = zeros(1,length(snr_dB), 'double');

% Needed later to record the frequency domain signal
video = VideoWriter('OFDM_Signal.avi');
video.FrameRate = 30;
open(video);

%Repitition coding reps
repFactor3 = 3;
repFactor5 = 5;

% Reading the video 
v = VideoReader('test_video.mp4');
frames = read(v,[1 1]); % reading the first 200 frames of the video

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


    % Measure bits in each frame
     num_bits = length(bits);

    %Repitition Channel Encoding
      txRepData3 = repelem(bits',repFactor3);
      txRepData5 = repelem(bits',repFactor5);

    %Perform 16-QAM modulation on every channel encoded bit stream
      
       bits_per_symb_qam = 4;

       %Group bits into 4
       grouped_qam = reshape(bits,bits_per_symb_qam,[]).';
    
       % Map each group of 4 bits to a 16-QAM symbol
       qam_mod = [];
       qam_mod = gray_map(bi2de(reshape(bits,4,[]).','left-msb')+1) - 7;
    
       qam_mod_3 = [];
       qam_mod_3 = gray_map(bi2de(reshape(txRepData3,4,[]).','left-msb')+1) - 7;
    
    
       qam_mod_5 = [];
       qam_mod_5 = gray_map(bi2de(reshape(txRepData5,4,[]).','left-msb')+1) - 7;

    % Convert serial data to parallel data

      parallel_data = reshape(qam_mod, [], 64);
      parallel_data_3 = reshape(qam_mod_3, [], 64);
      parallel_data_5 = reshape(qam_mod_5, [], 64);

     

        
    % Convert frequency domain signal to time domain
      ofdmSymbols_qam = ifft(parallel_data,[],2);
      ofdmSymbols_qam_3 = ifft(parallel_data_3,[],2);
      ofdmSymbols_qam_5 = ifft(parallel_data_5,[],2);



    % Add cyclic Prefix

    size_qam = size(ofdmSymbols_qam);
    ofdm_qam = zeros(size_qam(1),74);
    for num = 1:size(ofdmSymbols_qam, 1)
        data = ofdmSymbols_qam(num, :);
        cyclic_prefix = data (end-10+1:end);
        data_with_prefix= [cyclic_prefix data];
        ofdm_qam(num, :) = data_with_prefix;
    end


     size_qam_3 = size(ofdmSymbols_qam_3);
    ofdm_qam_3 = zeros(size_qam_3(1),74);
    for num = 1:size(ofdmSymbols_qam_3, 1)
        data_3 = ofdmSymbols_qam_3(num, :);
        cyclic_prefix_3 = data_3 (end-10+1:end);
        data_with_prefix_3= [cyclic_prefix_3 data_3];
        ofdm_qam_3(num, :) = data_with_prefix_3;
    end
    

      size_qam_5 = size(ofdmSymbols_qam_5);
    ofdm_qam_5 = zeros(size_qam_5(1),74);
    for num = 1:size(ofdmSymbols_qam_5, 1)
        data_5 = ofdmSymbols_qam_5(num, :);
        cyclic_prefix_5 = data_5 (end-10+1:end);
        data_with_prefix_5= [cyclic_prefix_5 data_5];
        ofdm_qam_5(num, :) = data_with_prefix_5;
    end
    
%% channel
    %Create Channels for every repetition coded stream

    numSymbols = length(ofdm_qam);
    QAM = unique(qam_mod);
    QAM_SNR_db = snr_dB(loop);
    QAM_symbols = QAM;
    % calculate the magnitude of each complex number
    QAM_magnitudeArray = abs(QAM_symbols);
    % square the magnitude of each complex number and sum them
    QAM_sumOfSquares = sum(QAM_magnitudeArray.^2);
    % divide the sum by the length of the original array
    QAM_c = QAM_sumOfSquares / length(QAM_symbols);
    QAM_SNR_W = db2pow(QAM_SNR_db);
    QAM_var_w = QAM_c/QAM_SNR_W;
    qam_txSig = ofdm_qam;
    taps=10;  %number of channel taps
    qam_channel = [randn(1,taps)+j*randn(1,taps) zeros(1,size(qam_txSig,2)-taps)].' ; % The Rayleigh fading channel padded with zeros
    qam_channel_mat = toeplitz(qam_channel); %  Convolution is the same as multiplying with the toeplitz matrix.
   
    
     QAM_3 = unique(qam_mod_3);
    QAM_SNR_db = snr_dB(loop);
    QAM_symbols_3 = QAM_3;
    % calculate the magnitude of each complex number
    QAM_magnitudeArray_3 = abs(QAM_symbols_3);
    % square the magnitude of each complex number and sum them
    QAM_sumOfSquares_3 = sum(QAM_magnitudeArray_3.^2);
    % divide the sum by the length of the original array
    QAM_c_3 = QAM_sumOfSquares_3 / length(QAM_symbols_3);
    QAM_SNR_W_3 = db2pow(QAM_SNR_db);
    QAM_var_w_3 = QAM_c_3/QAM_SNR_W_3;
    qam_txSig_3 = ofdm_qam_3;
    taps=10;  %number of channel taps
    qam_channel_3 = [randn(1,taps)+j*randn(1,taps) zeros(1,size(qam_txSig_3,2)-taps)].' ; % The Rayleigh fading channel padded with zeros
    qam_channel_mat_3 = toeplitz(qam_channel_3); %  Convolution is the same as multiplying with the toeplitz matrix.
    



     QAM_5 = unique(qam_mod_5);
    QAM_SNR_db = snr_dB(loop);
    QAM_symbols_5 = QAM_5;
    % calculate the magnitude of each complex number
    QAM_magnitudeArray_5 = abs(QAM_symbols_5);
    % square the magnitude of each complex number and sum them
    QAM_sumOfSquares_5 = sum(QAM_magnitudeArray_5.^2);
    % divide the sum by the length of the original array
    QAM_c_5 = QAM_sumOfSquares_5 / length(QAM_symbols_5);
    QAM_SNR_W_5 = db2pow(QAM_SNR_db);
    QAM_var_w_5 = QAM_c_5/QAM_SNR_W_5;
    qam_txSig_5 = ofdm_qam_5;
    taps=10;  %number of channel taps
    qam_channel_5 = [randn(1,taps)+j*randn(1,taps) zeros(1,size(qam_txSig_5,2)-taps)].' ; % The Rayleigh fading channel padded with zeros
    qam_channel_mat_5 = toeplitz(qam_channel_5); %  Convolution is the same as multiplying with the toeplitz matrix.
    




    % Pass all three signals through the channel

    qam_rcvsig=qam_channel_mat*qam_txSig'; % Pass through the Rayleigh fading channel
    qam_sig = qam_rcvsig + sqrt(QAM_var_w/2)*(randn(size(qam_rcvsig)) + 1i*randn(size(qam_rcvsig)));
    qam_sig_equal = inv(qam_channel_mat'*qam_channel_mat)*(qam_channel_mat')*qam_sig; % 0 forcing channel normalization
    y_qam = qam_sig_equal'; 


     qam_rcvsig_3=qam_channel_mat_3*qam_txSig_3'; % Pass through the Rayleigh fading channel
    qam_sig_3 = qam_rcvsig_3 + sqrt(QAM_var_w_3/2)*(randn(size(qam_rcvsig_3)) + 1i*randn(size(qam_rcvsig_3)));
    qam_sig_equal_3 = inv(qam_channel_mat_3'*qam_channel_mat_3)*(qam_channel_mat_3')*qam_sig_3; % 0 forcing channel normalization
    y_qam_3 = qam_sig_equal_3';


     qam_rcvsig_5=qam_channel_mat_5*qam_txSig_5'; % Pass through the Rayleigh fading channel
    qam_sig_5 = qam_rcvsig_5 + sqrt(QAM_var_w_5/2)*(randn(size(qam_rcvsig_5)) + 1i*randn(size(qam_rcvsig_5)));
    qam_sig_equal_5 = inv(qam_channel_mat_5'*qam_channel_mat_5)*(qam_channel_mat_5')*qam_sig_5; % 0 forcing channel normalization
    y_qam_5 = qam_sig_equal_5';
     


     

  %% Receiver
    

    % Remove cylcic prefix for rep = 1
    qam_received_cp_removed = y_qam(:,L+ 1:end);

    %Convert signal from time to frequency domain
    qam_received = fft(qam_received_cp_removed, [], 2);

    
    
   
       % Plot video for signal before and after normalization
         
 ofdm_qam_iv = ofdm_qam.';

for k = 1:numSymbols

    %Plot signal at OFDM Tx
    subplot(3, 1, 1);
    stem(1:74, abs(ofdm_qam_iv(:, k)), 'filled');
    xlabel('Subcarrier index');
    ylabel('Amplitude');
    title('OFDM Signal in Frequency Domain at Tx');

    % Plot received signal before equalization
    subplot(3, 1, 2);
    stem(1:74, abs(qam_sig(:, k)), 'filled');
    xlabel('Subcarrier index');
    ylabel('Amplitude');
    title('OFDM Signal in Frequency Domain (Before Equalization) ');

    % Plot received signal after equalization
    subplot(3, 1, 3);
    stem(1:74, abs(qam_sig_equal(:, k)), 'filled');
    xlabel('Subcarrier index');
    ylabel('Amplitude');
    title('OFDM Signal in Frequency Domain (After Equalization) ');

    frame = getframe(gcf);
    writeVideo(video, frame);
end 

close(video);

        

       

    %Convert parallel data to serial data
    serial_data = reshape(qam_received, [], 1);
    
    

    % Perform 16-QAM demodulation
     % 16 QAM
   QAM_bits_Rx = zeros(num_bits,1);
   I_bits = real(serial_data) > 0;
   Q_bits = imag(serial_data) > 0;
   I_mag = abs(real(serial_data));
   Q_mag = abs(imag(serial_data));
  
   QAM_bits_Rx(1:4:end) = I_bits;
   QAM_bits_Rx(2:4:end) = I_mag > 2/sqrt(10);
   QAM_bits_Rx(3:4:end) = Q_bits;
   QAM_bits_Rx(4:4:end) = Q_mag > 2/sqrt(10);
   qam_demod = QAM_bits_Rx';


   % Remove Cylcic Prefix for reps = 3
    qam_received_cp_removed_3 = y_qam_3(:,L+ 1:end);

    % Convert frequency domain signal to time domain
    qam_received_3 = fft(qam_received_cp_removed_3, [], 2);

    % Convert parallel data to serial data
    serial_data_3 = reshape(qam_received_3, [], 1);

    % Demodulate signal
   QAM_bits_Rx_3 = zeros(num_bits,1);
   I_bits_3 = real(serial_data_3) > 0;
   Q_bits_3 = imag(serial_data_3) > 0;
   I_mag_3 = abs(real(serial_data_3));
   Q_mag_3 = abs(imag(serial_data_3));
  
   QAM_bits_Rx_3(1:4:end) = I_bits_3;
   QAM_bits_Rx_3(2:4:end) = I_mag_3 > 2/sqrt(10);
   QAM_bits_Rx_3(3:4:end) = Q_bits_3;
   QAM_bits_Rx_3(4:4:end) = Q_mag_3 > 2/sqrt(10);
   qam_demod_3 = QAM_bits_Rx_3';


     % Remove Cylcic Prefix for reps = 5
    qam_received_cp_removed_5 = y_qam_5(:,L+ 1:end);

     % Convert frequency domain signal to time domain
    qam_received_5 = fft(qam_received_cp_removed_5, [], 2);

    % Convert parallel data to serial data
    serial_data_5 = reshape(qam_received_5, [], 1);


    % Demodulate signal
   QAM_bits_Rx_5 = zeros(num_bits,1);
   I_bits_5 = real(serial_data_5) > 0;
   Q_bits_5 = imag(serial_data_5) > 0;
   I_mag_5 = abs(real(serial_data_5));
   Q_mag_5 = abs(imag(serial_data_5));
  
   QAM_bits_Rx_5(1:4:end) = I_bits_5;
   QAM_bits_Rx_5(2:4:end) = I_mag_5 > 2/sqrt(10);
   QAM_bits_Rx_5(3:4:end) = Q_bits_5;
   QAM_bits_Rx_5(4:4:end) = Q_mag_5 > 2/sqrt(10);
   qam_demod_5 = QAM_bits_Rx_5';
    

     % Perform Channel Decoding
        blocks3 = reshape(qam_demod_3, repFactor3, []).';
        blocks5 = reshape(qam_demod_5, repFactor5, []).';


        ones_count3 = sum(blocks3, 2);
        ones_count5 = sum(blocks5, 2);
        
        % Decode the bit stream
        decoded_bits3 = ones_count3 >= repFactor3/2;
        decoded_bits3 = decoded_bits3(:)';

        decoded_bits5 = ones_count5 >= repFactor5/2;
        decoded_bits5 = decoded_bits5(:)';



    
    % Recreate the video for rep = 1
         orig_class = class(currentFrame);
     orig_size = size(currentFrame);
     qam_reconstructedFrames(:,:,:,i) = reshape(typecast(uint8(bin2dec ...
         (char(reshape(qam_demod, 8, [])+'0').')),orig_class), orig_size);


    % Recreate the video for rep = 3
     orig_class = class(currentFrame);
     orig_size = size(currentFrame);
     qam_reconstructedFrames_3(:,:,:,i) = reshape(typecast(uint8(bin2dec ...
         (char(reshape(decoded_bits3, 8, [])+'0').')),orig_class), orig_size);


      % Recreate the video for rep = 5
          orig_class = class(currentFrame);
     orig_size = size(currentFrame);
     qam_reconstructedFrames_5(:,:,:,i) = reshape(typecast(uint8(bin2dec ...
         (char(reshape(decoded_bits5, 8, [])+'0').')),orig_class), orig_size);


% Calculate BER 
    BER_qam(loop) = BER_qam(loop) + sum(qam_demod ~= bits.');
    BER_qam_3(loop) = BER_qam_3(loop) + sum(decoded_bits3 ~= bits);
    BER_qam_5(loop) = BER_qam_5(loop) + sum(decoded_bits5 ~= bits);
    %disp(i)

end

% Normalize BERs over frames
BER_qam(loop) = BER_qam(loop)/(num_bits*size_frame);
BER_qam_3(loop) = BER_qam_3(loop)/(num_bits*size_frame);
BER_qam_5(loop) = BER_qam_5(loop)/(num_bits*size_frame);


disp(loop)
end 


% Plot BER vs. SNR for rep = 1
figure(1)
semilogy(snr_dB, BER_qam);
xlabel('SNR (dB)')
ylabel('16-QAM BER (Pe)')
title('16-QAM Probability of error(1 repitition) vs. SNR')


% Plot BER vs. SNR for rep = 3
figure(2)
semilogy(snr_dB, BER_qam_3);
xlabel('SNR (dB)')
ylabel('16-QAM BER (Pe)')
title('16-QAM Probability of error(3 repitition) vs. SNR')


% Plot BER vs. SNR for rep = 5
figure(3)
semilogy(snr_dB, BER_qam_5);
xlabel('SNR (dB)')
ylabel('16-QAM BER (Pe)')
title('16-QAM Probability of error(5 repitition) vs. SNR')

% Play the video
implay(qam_reconstructedFrames);