%% Convolutional Code

Rc=0.5; % Code rate 
const_length=3; % Constraint length of the convolutional code
g1=7; % First octal generator 
g2=5; % Second octal generator
tb=const_length-1; % Terminating bits added to end of message (equal to number of memory elements) to force the encoder to contain all zero bits in the memory elements
K=1000-tb; % Message length (the added terminating bits do not contain any information)
N=(K+tb)/Rc; % Codeword length with terminating bits 
message=zeros(1,K+tb); 
codeword=zeros(1,N+(tb/Rc)); % Number of terminating coded bits added to end of codeword 
snr=[0 1 2 3 4 5 6 7 8 9 10];
frames=[10 100 100 100 100 1000 5000 10000 10000 10000 10000]; % Number of codewords generated
snr_uncoded=linspace(0,max(snr),length(snr));
Eb_N0=10.^(snr/10);
ber=zeros(1,length(snr));
trellis=poly2trellis(const_length,[g1 g2]); % Matlab constructs a trellis based on the octal generators and constraint length of the convolutional code

for s=1:length(snr)
    errors=0; % Error counter
    sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
    for f=1:frames(s)
        message=randi([0 1],1,K); % Generate random 0 and 1
                
        for i=1:tb  
            message(K+i)=0; % Add terminating bits 
        end 
        
        codeword=convenc(message,trellis); % Convolutional encoder  
        x=1-2*codeword;
        n=sigma*randn(1,length(x)); % Generate AWGN 
        y=x+n; % Add noise to BPSK signal
        d=(y<0);
        decoded_message=vitdec(d,trellis,5*const_length,'term','hard'); % Hard-decision Viterbi decoding 
         
        for i=1:K  
            if decoded_message(i)~=message(i) 
                errors=errors+1; % Count errors
            end
        end
    end
    ber(s)=errors/(frames(s)*K); % Calculate BER
end

ber_uncoded=0.5*erfc(sqrt(Eb_N0));

% Enhanced visualization
figure(1);
semilogy(snr_uncoded, ber_uncoded, '-o', 'LineWidth', 1, 'MarkerSize', 8); hold on;
semilogy(snr, ber, '-s', 'LineWidth', 1, 'MarkerSize', 8);
grid on;

% Annotations
legend1 = 'Uncoded BPSK modulation';
legend2 = sprintf('(%d,%d) convolutional code',g1,g2);
legend({legend1, legend2}, 'Location', 'SouthWest', 'FontSize', 10);
title('BER Performance of Convolutional Code and Uncoded BPSK', 'FontSize', 12);
xlabel('Eb/N0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);

hold off;

%% Convolutional Codes with different Polynomials

% Common Parameters
Rc=0.5; % Code rate
snr=[0 1 2 3 4 5 6 7 8 9 10]; % SNR in dB
Eb_N0=10.^(snr/10); % Convert SNR to linear scale
frames=[10 100 100 100 100 1000 5000 10000 10000 10000 10000]; % Number of codewords generated
ber_uncoded = 0.5 * erfc(sqrt(Eb_N0)); % Uncoded BPSK BER

% Define convolutional code configurations
convolutional_configs = {
    struct('const_length', 3, 'generators', [7,5], 'label', '(7,5)_8 Convolutional Code');
    struct('const_length', 4, 'generators', [17,15], 'label', '(17,15)_8 Convolutional Code');
    struct('const_length', 7, 'generators', [133,171], 'label', '(133,171)_8 Convolutional Code'),
};

% Plot uncoded BPSK for comparison
figure;
semilogy(snr, ber_uncoded, '-o','LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Uncoded BPSK') % Plot uncoded BPSK for comparison
hold on;

% Loop through each convolutional code configuration
for config_i = 1:length(convolutional_configs)
    config = convolutional_configs{config_i};
    const_length = config.const_length;
    generators = config.generators;
    trellis = poly2trellis(const_length, generators); % Generate a trellis
    tb = const_length - 1; % Terminating bits
    K = 1000 - tb; % Message length
    ber = zeros(1, length(snr)); % Initializing BER array

    % Simulate BER for each SNR
    for s = 1:length(snr)
        errors=0; % Error counter
        sigma=sqrt(1/(2*Rc*Eb_N0(s))); % Standard deviation of noise
        for f=1:frames(s)
            message=randi([0 1],1,K); % Generate random 0 and 1
            for i=1:tb  
                message(K+i)=0; % Add terminating bits 
            end
            codeword=convenc(message,trellis); % Convolutional encoder  
            x=1-2*codeword; % BPSK modulated signal
            n=sigma*randn(1,length(x)); % Generate AWGN 
            y=x+n; % Add noise to BPSK signal
            d=(y<0); % Demodulate signal
            decoded_message=vitdec(d, trellis, 5*const_length,'term','hard'); % Hard-decision Viterbi decoding
            for i=1:K  
                if decoded_message(i)~=message(i) 
                    errors=errors+1; % Count errors
                end
            end
        end
        ber(s)=errors/(frames(s)*K); % Calculate BER
    end

    % Plot BER for this configuration
    semilogy(snr, ber, '-o', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', config.label);
end

% Final Plot
legend('show', 'Location', 'southwest', 'FontSize', 10);
xlabel('Eb/N0, dB', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER Performance of Different Convolutional Codes', 'FontSize', 12);
grid on;
hold off;