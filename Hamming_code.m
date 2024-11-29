%% Hamming Code

K=4; % Message length 
N=7; % Codeword length 
D=3; % Minimum Hamming distance
Rc=K/N; % Code rate 
message=zeros(1,K);
codeword=zeros(1,N); 
snr=[0 1 2 3 4 5 6 7 8 9 10];
Eb_N0=10.^(snr/10);
frames=[10000 10000 10000 10000 10000 10000 100000 100000 100000 1000000 1000000]; % Number of codewords generated
snr_uncoded=linspace(0,max(snr),length(snr));
ber=zeros(1,length(snr));
G=[[1,0,0,0,1,1,0];[0,1,0,0,0,1,1];[0,0,1,0,1,1,1];[0,0,0,1,1,0,1]]; % Define the generator matrix of the (7, 4, 3) Hamming code
H=[[1,0,1,1,1,0,0];[1,1,1,0,0,1,0];[0,1,1,1,0,0,1]]; % Define the parity-check matrix of the (7, 4, 3) Hamming code
synd_table=[[1 1 0];[0 1 1];[1 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1]]; % Define the syndrome table for each possible 1 bit-error in the 7-bit codeword

for s=1:length(snr)
    errors=0; % Error counter
   
    sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
    for f=1:frames(s)
        message=randi([0 1],1,K); % Generate random k-bit message  
        codeword=mod(message*G,2); % Encode using generator matrix
        x=1-2*codeword; % BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
    
        n=sigma*randn(1,length(x)); % Generate AWGN 
        y=x+n; % Add noise to BPSK signal
     
        r=(y<0); % BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0
     
        error_pattern=zeros(1,N); % Initialise error pattern vector to all zeros
        syndromes=mod(H*r',2); % Calculate the syndromes by multiplying the parity-check matrix by the transpose of the received binary word
        

        for i=1:N
            a=isequal(syndromes',synd_table(i,:)); % Compare syndrome vector with every possible syndrome in the syndrome table
            if a==1 % If the syndrome matches the i-th syndrome in the syndrome table
                error_pattern(i)=1; % The bit error it in the i-th bit of the received binary word r
            end
        end
        decoded_message=mod(r+error_pattern,2); % Add the error pattern to the received binary word r

        for i=1:K  
            if decoded_message(i)~=message(i) 
                errors=errors+1; % Count errors
            end
        end
     end
     ber(s)=errors/(frames(s)*K); % Calculate BER
end

ber_uncoded=0.5*erfc(sqrt(Eb_N0)); % The BER for BPSK

% Enhanced visualization
figure(1);
semilogy(snr_uncoded, ber_uncoded, '-o', 'LineWidth', 1, 'MarkerSize', 8); hold on;
semilogy(snr, ber, '-s', 'LineWidth', 1, 'MarkerSize', 8);
grid on;

% Annotations
legend1 = 'Uncoded BPSK modulation';
legend2 = sprintf('(%d,%d,%d) Hamming code',N,K,D);
legend({legend1, legend2}, 'Location', 'SouthWest', 'FontSize', 10);
title('BER Performance of Hamming Code and Uncoded BPSK', 'FontSize', 12);
xlabel('Eb/N0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);

hold off;

%% Hamming Codes Vs Repetition Codes

%Common parameters
snr=[0 1 2 3 4 5 6 7 8 9 10];
Eb_N0=10.^(snr/10); %Convert SNR to linear scale
frames=[10000 10000 10000 10000 10000 10000 100000 100000 100000 1000000 1000000]; %number of codewords generated
ber_uncoded = 0.5 * erfc(sqrt(Eb_N0)); % BER for uncoded BPSK

% List of repetition Lengths
repetition_lengths = [3,5,101,1001];

% BER storage for repetition codes
ber_repetition = zeros(length(repetition_lengths), length(snr));

% Plot uncoded BPSK for comparision
figure;
semilogy(snr, ber_uncoded, '-o','LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Uncoded BPSK') % Plot uncoded BPSK for comparison
hold on;

% Repetition code Simulation
for rep = 1:length(repetition_lengths)
    N = repetition_lengths(rep); % Codeword length
    K = 1; % Message Length
    D = N; % Minimum Hamming Distance
    Rc = K / N; % Code Rate
    ber = zeros(1, length(snr)); % Initialize BER array
    G = ones(1, N); % Generator Matrix for Repetition code

    % BER simulation for repetition codes
    for s = 1:length(snr)
        errors = 0; % Error Counter
        sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
        for f = 1:frames(s)
            message = randi([0 1], 1, K); % Random message
            codeword = mod(message * G,2); % Encode message using Generator Matrix
            x=1-2*codeword; % BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
            n=sigma*randn(1,length(x)); % Generate AWGN 
            y=x+n; % Add noise to BPSK signal
            r=(y<0); % BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0
            n1 = nnz(r); % Count 1s in received codeword
            n0 = N - n1; % Count 0s in received codeword
            decoded_message = n1 > n0; % Majority Decoding
            if decoded_message ~= message
                errors = errors + 1;
            end
        end
        ber(s) = errors / (frames(s) * K); % Calculate BER
    end
    ber_repetition(rep, :) = ber; % Store BER results
    semilogy(snr, ber, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', sprintf('(%d,%d,%d) Repetition Code', N, K, D));
end

% Hamming code Simulation
% Parameters for Hamming Code
K = 4; % Message length
N = 7; % Codeword length
D = 3; % Minimum Hamming distane
Rc = K / N; % Code rate
G=[[1,0,0,0,1,1,0];[0,1,0,0,0,1,1];[0,0,1,0,1,1,1];[0,0,0,1,1,0,1]]; % Define the generator matrix of the (7, 4, 3) Hamming code
H=[[1,0,1,1,1,0,0];[1,1,1,0,0,1,0];[0,1,1,1,0,0,1]]; % Define the parity-check matrix of the (7, 4, 3) Hamming code
synd_table=[[1 1 0];[0 1 1];[1 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1]]; % Define the syndrome table for each possible 1 bit-error in the 7-bit codeword
ber_hamming = zeros(1, length(snr)); % Initialize BER array

% BER Simulation for Hamming Code
for s=1:length(snr)
    errors=0; % Error counter
    sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
    for f=1:frames(s)
        message = randi([0 1],1,K); % Generate random K-bit message  
        codeword = mod(message*G,2); % Encode message using generator matrix
        x=1-2*codeword; % BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
        n=sigma*randn(1,length(x)); % Generate AWGN 
        y=x+n; % Add noise to BPSK signal
        r=(y<0); % BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0
        error_pattern = zeros(1,N); % Initialise error pattern vector to all zeros
        syndromes = mod(H*r',2); % Calculate the syndromes by multiplying the parity-check matrix by the transpose of the received binary word
        for i=1:N
            a=isequal(syndromes',synd_table(i,:)); % Compare syndrome vector with every possible syndrome in the syndrome table
            if a==1 % If the syndrome matches the i-th syndrome in the syndrome table
                error_pattern(i) = 1; % The bit error it in the i-th bit of the received binary word r
            end
        end
        decoded_message = mod(r+error_pattern,2); % Add the error pattern to the received binary word r
        for i=1:K  
            if decoded_message(i)~=message(i) 
                errors=errors+1; % Count errors
            end
        end
     end
     ber_hamming(s) = errors/(frames(s)*K); % Calculate BER
end

% Plot Hamming Code BER
semilogy(snr, ber_hamming, '-s', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', sprintf('(%d,%d,%d) Hamming Code', N, K, D));

% Finalize Plot
legend('show', 'Location', 'southwest', 'FontSize', 10);
xlabel('Eb/N0, dB', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER Performance of Hamming Code and Repetition Codes', 'FontSize', 12);
grid on;
hold off;