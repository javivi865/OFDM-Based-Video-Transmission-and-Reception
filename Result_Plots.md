# Results and Plots

## Probability of Error vs. Received SNR
![BER_vs_SNR](https://github.com/javivi865/OFDM-Based-Video-Transmission-and-Reception/assets/169334804/a912bc60-2ae5-477c-8574-98268215a81a)

## Probability of Error vs. Received SNR for 16-QAM with Repetition Coding
![BER_vs_SNR_16QAM_Repetiton_Coding](https://github.com/javivi865/OFDM-Based-Video-Transmission-and-Reception/assets/169334804/a2fc56b1-dde4-48f6-8e5e-34c9ea58aeec)

Observation: The bit error rate reduces with the number of repetitions using the repetition coding technique. For a 16-QAM system, for a single repetition at 10dB SNR, we get a BER of 0.25; for 3 repetitions of each bit, BER is 0.19, and for 5 repetitions, BER is 0.14. Channel encoding reduces the probability of error for a given SNR.

## Bit Rate Calculations
If a given Bandwidth is given, the bit rate can be calculated easily using the following formula: Bit rate = Bandwidth * log2(N). If the system has a 10MHz bandwidth, the bit rate for the BPSK, QPSK, and 16-QAM are:

For BPSK, Bit rate = 10e6 * log2(1) = 10Mbps

For QPSK, Bit rate = 10e6 * log2(4) = 20Mbps

For 16-QAM, Bit rate = 10e6 * log2(16) = 40Mbps

