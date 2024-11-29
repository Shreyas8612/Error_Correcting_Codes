%% Soft Decision Decoding

% Common Parameters
snr=[0 1 2 3 4 5 6 7 8 9 10]; % SNR in dB
Eb_N0=10.^(snr/10); % Convert SNR to linear scale
ber_uncoded = 0.5 * erfc(sqrt(Eb_N0)); % Uncoded BPSK BER

% Plot uncoded BPSK for comparison
figure;
semilogy(snr, ber_uncoded, '-s','LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Uncoded BPSK') % Plot uncoded BPSK for comparison
hold on;

%% Repetition Code (3,1,3)
% Parameters for repetition code
N_rep = 3; % Codeword length
K_rep = 1; % Message length
frames_rep=[5000 5000 5000 5000 5000 5000 50000 50000 50000 50000 50000]; % Number of codewords generated
Rc_rep = K_rep / N_rep; % Code rate
G_rep = ones(1, N_rep); % Generator matrix
ber_rep = zeros(1, length(snr)); % Initialize BER

% Simulate Repetition code
for s = 1:length(snr)
    errors = 0; % Initialize error counter
    sigma = sqrt(1/(2 * Rc_rep * Eb_N0(s))); % Standard deviation of noise
    
    % Loop over frames for each SNR
    for f = 1:frames_rep(s)
        message = randi([0 1], 1, K_rep); % Generate 1-bit random message
        codeword = mod(message * G_rep, 2); % Repetition encoding
        x = 1 - 2 * codeword; % BPSK modulation
        n = sigma * randn(1, length(x)); % Gaussian noise
        y = x + n; % Add noise to BPSK signal
        
        % Decoding by majority
        d = (y < 0); % Demodulation: positive -> 0, negative -> 1
        n1 = nnz(d); % Count 1s in received codeword
        n0 = N_rep - n1; % Count 0s in received codeword
        decoded_message = n1 > n0; % Majority decoding
        if decoded_message ~= message % If decoded incorrectly
            errors = errors + 1;
        end
    end
    ber_rep(s) = errors / (frames_rep(s) * K_rep); % Calculate BER for SNR
end

% Plot repetition code BER
semilogy(snr, ber_rep, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', '(3,1,3) Repetition Code');

%% Hamming Code (7,4,3)
% Parameters for Hamming Code
K_ham = 4; % Message length
N_ham = 7; % Codeword length
D_ham = 3; % Minimum Hamming distane
Rc_ham = K_ham / N_ham; % Code rate
frames_ham=[2000 2000 2000 2000 2000 2000 20000 20000 20000 200000 2000000]; % Number of codewords generated
G_ham=[[1,0,0,0,1,1,0];[0,1,0,0,0,1,1];[0,0,1,0,1,1,1];[0,0,0,1,1,0,1]]; % Define the generator matrix of the (7, 4, 3) Hamming code
H_ham=[[1,0,1,1,1,0,0];[1,1,1,0,0,1,0];[0,1,1,1,0,0,1]]; % Define the parity-check matrix of the (7, 4, 3) Hamming code
synd_table_ham=[[1 1 0];[0 1 1];[1 1 1];[1 0 1];[1 0 0];[0 1 0];[0 0 1]]; % Define the syndrome table for each possible 1 bit-error in the 7-bit codeword
ber_ham = zeros(1, length(snr)); % Initialize BER array

% BER Simulation for Hamming Code
for s=1:length(snr)
    errors=0; % Error counter
    sigma=sqrt(1/(2*Rc_ham*Eb_N0(s))); % Standard deviation of noise
    for f=1:frames_ham(s)
        message = randi([0 1],1,K_ham); % Generate random K-bit message  
        codeword = mod(message*G_ham,2); % Encode message using generator matrix
        x=1-2*codeword; % BPSK signal obtained by mapping a bit 0 to +1 (in-phase carrier signal) and a bit 1 to -1 (out-of-phase carrier signal)
        n=sigma*randn(1,length(x)); % Generate AWGN 
        y=x+n; % Add noise to BPSK signal
        d=(y<0); % BPSK demodulation. If the received signal is negative then the demodulated bit is 1, else the demodulated bit is 0
        error_pattern = zeros(1,N_ham); % Initialise error pattern vector to all zeros
        syndromes = mod(H_ham*d',2); % Calculate the syndromes by multiplying the parity-check matrix by the transpose of the received binary word
        for i=1:N_ham
            a=isequal(syndromes',synd_table_ham(i,:)); % Compare syndrome vector with every possible syndrome in the syndrome table
            if a==1 % If the syndrome matches the i-th syndrome in the syndrome table
                error_pattern(i) = 1; % The bit error it in the i-th bit of the received binary word r
            end
        end
        decoded_message = mod(d+error_pattern,2); % Add the error pattern to the received binary word r
        for i=1:K_ham  
            if decoded_message(i)~=message(i) 
                errors=errors+1; % Count errors
            end
        end
     end
     ber_ham(s) = errors/(frames_ham(s)*K_ham); % Calculate BER
end

% Plot Hamming Code BER
semilogy(snr, ber_ham, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', '(7,4,3) Hamming Code');

%% Convolutional Codes Comparison
% Define convolutional code configurations
convolutional_configs = {
    struct('const_length', 3, 'generators', [7, 5], 'label', '(7,5)_8 Convolutional Code')...
    struct('const_length', 4, 'generators', [17, 15], 'label', '(17,15)_8 Convolutional Code', 'soft', true)...
    struct('const_length', 7, 'generators', [133, 171], 'label', '(133,171)_8 Convolutional Code'),
};

Rc_conv = 0.5; % Code rate
frames_conv = [1000 1000 1000 1000 1000 5000 10000 20000 20000 50000 50000]; % Number of codewords generated

% Loop through each convolutional code configuration
for config_i = 1:length(convolutional_configs)
    config = convolutional_configs{config_i};
    const_length = config.const_length;
    generators = config.generators;
    trellis = poly2trellis(const_length, generators); % Generate a trellis
    tb = const_length - 1; % Terminating bits
    K_conv = 1000 - tb; % Message length
    ber_conv = zeros(1, length(snr)); % Initializing BER array for hard decision
    if isfield(config, 'soft') && config.soft
        ber_conv_soft = zeros(1, length(snr)); % Initializing BER array for soft decision
    end

    % Simulate BER for each SNR
    for s = 1:length(snr)
        errors=0; % Error counter
        errors_soft = 0; % Error counter for soft decision
        sigma = sqrt(1 / (2 * Rc_conv * Eb_N0(s))); % Standard deviation of noise
        for f = 1:frames_conv(s)
            message = randi([0 1], 1, K_conv); % Generate random 0 and 1
            for i=1:tb  
                message(K_conv + i) = 0; % Add terminating bits 
            end
            codeword = convenc(message,trellis); % Convolutional encoder  
            x = 1 - 2 * codeword; % BPSK modulated signal
            n = sigma * randn(1, length(x)); % Generate AWGN 
            y = x + n; % Add noise to BPSK signal

            % Hard-decision Viterbi decoding
            d = (y < 0); % Demodulate signal
            decoded_hard = vitdec(d, trellis, 5*const_length,'term','hard');
            for i=1:K_conv  
                if decoded_hard(i) ~= message(i) 
                    errors = errors + 1; % Count errors
                end
            end

            % Soft-decision Viterbi decoding
            if isfield(config, 'soft') && config.soft
                decoded_soft = vitdec(y, trellis, 5 * const_length, 'term', 'unquant');
                for i=1:K_conv
                    if decoded_soft(i) ~= message(i) 
                    errors_soft = errors_soft + 1; % Count errors
                    end
                end
            end
        end
        ber_conv(s) = errors / (frames_conv(s) * K_conv); % Calculate BER
        if isfield(config, 'soft') && config.soft
            ber_conv_soft(s) = errors_soft / (frames_conv(s) * K_conv);
        end
    end

    % Plot BER for this configuration
    semilogy(snr, ber_conv, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', config.label);
    if isfield(config, 'soft') && config.soft
        semilogy(snr, ber_conv_soft, '-s', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', [config.label ' (Soft Decision']);
    end
end

%% Final Plot
legend('show', 'Location', 'southwest', 'FontSize', 10);
xlabel('Eb/N0, dB', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER Performance Comparison of Repetition, Hamming, and Convolutional Codes', 'FontSize', 12);
grid on;
hold off;
