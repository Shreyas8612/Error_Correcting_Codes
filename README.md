# Error-Correcting Codes: MATLAB Analysis

This repository contains MATLAB scripts and results for the analysis and comparison of different error-correcting codes: Repetition Codes, Hamming Codes, and Convolutional Codes. The work is aimed at evaluating their Bit-Error Rate (BER) performance under different conditions and comparing them to uncoded Binary Phase Shift Keying (BPSK).

## Overview
The MATLAB code in this repository is used to simulate and analyze the BER performance of different error-correcting codes for various Signal-to-Noise Ratios (SNR). Specifically, the comparison includes:

- **Repetition Code (3,1,3), (5,1,5), (101,1,101), (1001,1,1001)**
- **Hamming Code (7,4,3)**
- **Convolutional Codes** with different configurations (e.g., (7,5), (17,15), (133,171))

The BER performance is plotted against uncoded BPSK to observe the effectiveness of each error-correcting scheme.

## Files in This Repository
- **`Error_Correcting_Codes.m`**: Different MATLAB scripts to simulate and analyze the BER performance of Repetition Codes, Hamming Codes, and Convolutional Codes.

## Coding Schemes and Analysis
1. **Repetition Codes**
   - Uses redundancy by repeating each bit multiple times to improve error correction.
   - Effective for low noise scenarios but has reduced energy efficiency due to spreading the available energy across multiple repeated bits.
   - BER results show that the Repetition Codes have worse performance compared to uncoded BPSK at higher SNR due to reduced code rate.

2. **Hamming Code (7,4,3)**
   - Uses structured parity-check bits to enable single-bit error correction.
   - Initially shows worse BER at very low SNR, but outperforms uncoded BPSK at moderate to high SNR due to its error correction capabilities.
   - A good trade-off between redundancy and error correction efficiency.

3. **Convolutional Codes**
   - Simulated with different constraint lengths and generator polynomials, e.g., (7,5), (17,15), (133,171).
   - Performance improves with longer constraint lengths, resulting in more memory elements being used for encoding and better BER performance.
   - Soft decision decoding further improves the BER by considering the likelihood of received bits.

## Key Results
- **BER Performance**: The BER performance of each coding scheme was evaluated against uncoded BPSK. Repetition Codes performed the worst due to low energy efficiency, while Hamming Codes and Convolutional Codes showed significant improvement at higher SNR.
- **Soft vs Hard Decoding**: Convolutional Codes with soft decision decoding outperformed those with hard decision decoding, highlighting the benefits of analyzing the likelihood of received bits for better error correction.

## How to Run the Code
1. Clone the repository to your local machine:
   ```bash
   git clone git@github.com:Shreyas8612/Error_Correcting_Codes.git
   ```
2. Open MATLAB and navigate to the cloned directory.
3. Run the script `Error_Correcting_Codes.m` to simulate the BER performance of each code.
4. The results will be displayed in MATLAB figures, showcasing the BER vs. SNR performance for each scheme.

## Conclusion
This project provides an insightful comparison of different error-correcting codes, demonstrating their strengths and limitations in improving BER under various SNR conditions. Convolutional Codes, especially with soft decision decoding, showed the best performance, indicating their superior error correction capabilities compared to Repetition and Hamming Codes.

## Author
**Shreyas Ravi**

## License
This project is licensed under the MIT License.

