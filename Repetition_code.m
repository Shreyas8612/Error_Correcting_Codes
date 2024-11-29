%% Repetition Code

K=1; % Message length 
N=3; % Codeword length 
D=N; % Minimum Hamming distance
Rc=K/N; % Code rate
message=zeros(1,K); % Initialise k-bit message to all zeros 
codeword=zeros(1,N); % Initialise n-bit codeword to all zeros
snr=[0 1 2 3 4 5 6 7 8 9 10]; % Range of signal-to-noise ratios in dB
frames=[1000 1000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000]; % Number of codewords generated for each SNR
snr_uncoded=linspace(0,max(snr),length(snr)); % Range of signal-to-noise ratios for BER of uncoded BPSK which we will compare our results with
Eb_N0=10.^(snr/10); % Convert snr in dB to Eb/N0
ber=zeros(1,length(snr)); % Initialise vector of bit-error rates for each snr value
G=ones(1,N); % Repetition code generator matrix


for s=1:length(snr)
    errors=0; % Initialise error counter at the start of each snr
    sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
    for f=1:frames(s)
        message=randi([0 1],1,K); % Generate 1-bit message with random 0s and 1s    
        codeword=mod(message*G,2); % Multiply message by generator matrix to obtain the n-bit repetition codeword
        x=1-2*codeword; % BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
    
        n=sigma*randn(1,length(x)); % Generate white Gaussian noise
        y=x+n; % Add noise to BPSK signal
     
        d=(y<0); % BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0

        % Majority Decoding

        n1=nnz(d); % Count the number of 1s in the received binary word d
        n0=N-n1; % The number of 0s in the recevied binary word d

        if n1>n0 % If the number of 1s is greater than the number 0s in d
             decoded_message=1; % Decoded message bit is 1
        else
            decoded_message=0;  % Else the decoded message bit is 0
        end
         
        if decoded_message~=message % If the decoded message bit is different to the message bit then this is a bit error
                errors=errors+1; % Add the bit error to the error counter
        end
       
    end
    ber(s)=errors/(frames(s)*K); % Calculate BER of the repetition code
end

ber_uncoded=0.5*erfc(sqrt(Eb_N0)); % The BER for BPSK

% Enhanced visualization
figure(1);
semilogy(snr_uncoded, ber_uncoded, '-o', 'LineWidth', 1, 'MarkerSize', 8); hold on;
semilogy(snr, ber, '-s', 'LineWidth', 1, 'MarkerSize', 8);
grid on;

% Annotations
legend1 = 'Uncoded BPSK modulation';
legend2 = sprintf('(%d,%d,%d) Repetition Code', N, K, D);
legend({legend1, legend2}, 'Location', 'SouthWest', 'FontSize', 10);
title('BER Performance of Repetition Code and Uncoded BPSK', 'FontSize', 12);
xlabel('Eb/N0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);

hold off;

%% Repetition Codes with different lengths

K = 1; % Message length
snr = [0 1 2 3 4 5 6 7 8 9 10]; % Range of SNR in dB
frames = [1000 10000 10000 100000 1000000 1000000 1000000 1000000 1000000 1000000 1000000];
Eb_N0 = 10.^(snr/10); % Convert SNR to linear scale for calculation
ber_uncoded = 0.5 * erfc(sqrt(Eb_N0)); % BER for uncoded BPSK

% List of repetition lengths to test
repetition_lengths = [3, 5, 101, 1001];

figure;
semilogy(snr, ber_uncoded, '-o','LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Uncoded BPSK'); % Plot uncoded BPSK for comparison
hold on;

% Loop over different repetition lengths
for N = repetition_lengths
    D = N; % Minimum Hamming distance
    Rc = K / N; % Code rate
    ber = zeros(1, length(snr)); % Initialize BER for each SNR value
    G = ones(1, N); % Repetition code generator matrix
    
    % Loop over SNR values
    for s = 1:length(snr)
        errors = 0; % Initialize error counter
        sigma = sqrt(1/(2 * Rc * Eb_N0(s))); % Standard deviation of noise
        
        % Loop over frames for each SNR
        for f = 1:frames(s)
            message = randi([0 1], 1, K); % Generate 1-bit random message
            codeword = mod(message * G, 2); % Repetition encoding
            x = 1 - 2 * codeword; % BPSK modulation
            n = sigma * randn(1, length(x)); % Gaussian noise
            y = x + n; % Add noise to BPSK signal
            
            % Decoding by majority
            d = (y < 0); % Demodulation: positive -> 0, negative -> 1
            n1 = nnz(d); % Count 1s in received codeword
            n0 = N - n1; % Count 0s in received codeword
            decoded_message = n1 > n0; % Majority decoding
            if decoded_message ~= message % If decoded incorrectly
                errors = errors + 1;
            end
        end
        ber(s) = errors / (frames(s) * K); % Calculate BER for SNR
    end
    
    % Plot BER curve for this repetition length
    semilogy(snr, ber, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', sprintf('(%d,%d,%d) Repetition Code', N, K, D));
end

hold off;
legend('Location', 'southwest', 'FontSize', 10);
xlabel('Eb/N0, dB', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER Performance of Repetition Codes with Different Lengths', 'FontSize', 12);
grid on;