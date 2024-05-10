# OFDM-Based Video Transmission and Reception in MATLAB. 

This communication system includes one transmitter and one receiver with a single antenna. The transmitter supports BPSK, QPSK, and 16-QAM modulation. The OFDM system has 64 subcarriers. The channel between the transmitter and receiver has 10 taps; the channel at each tap follows a Rayleigh distribution. i.e., ℎ[𝑙]~ 𝐶𝑁(0,1). The channel is assumed to be time-invariant within each OFDM symbol. Further, the channels at the different OFDM symbols are assumed to be IID.

The communication system will be tested by sending video frames through the system. There will be 3 .m files corresponding to each modulation (BPSk, QPSk, 16-QAM). For all 3 modulations, there will be code to plot the Prob. of Error (BER) vs. SNR ranging from -10 dB to 10 dB. For the 16-QAM file, repetition coding with 1, 3, and 5 was added to test the system and attempt to improve


