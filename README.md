# OFDM-Based Video Transmission and Reception in MATLAB. 

### This communication system includes:
One transmitter and one receiver with a single antenna. The transmitter supports BPSK, QPSK, and 16-QAM modulation. The OFDM system has 64 subcarriers. The channel between the transmitter and receiver has 10 taps. The channel at each tap follows a Rayleigh distribution. i.e., ℎ[𝑙]~ 𝐶𝑁(0,1). The channel is assumed to be time-invariant within each OFDM symbol. Further, the channels at the different OFDM symbols are assumed to be IID.

### Testing the communication system: 
The communication system will be tested by sending video frames through the system. 3 MATLAB files (.m extension) will correspond to each modulation (BPSk, QPSk, 16-QAM). For all 3 modulations, there will be code to plot the Probability of Error (BER) vs. SNR ranging from -10 dB to 10 dB. For the 16-QAM file, repetition coding with 1, 3, and 5 to notice the different across this technique. In addition, the 16-QAM will display the OFDM signal in the frequency domain before and after the channel equalization. Lastly, every modulation will play the received frames of the video. 

### Repository Files
1- Video used for testing of the communication system. 

2- Results and Plots file.

3- MATLAB files for each modulation. 


