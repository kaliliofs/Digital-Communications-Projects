%% Matlab Intialization
clc;
clear;                                      % Clear matlab work screen

%% Signal Intializatin
n_sym = 1e6;                                % Number of symbols
s = randi([0 1], n_sym, 1);                 % Random bit generation

%% Encoding Intialization
% 0: Without encoding, 1: Repetition, 2: Linear Block, 3: Convlution
% 4: Changing constraint length in convlution code
operation = 0;                          % Select encoding method              

% Repetition code
rep_num = 11;                           % Number of repetition in CC

% Linear block code 
n = 7;                                  % Code length at a time
k = 3;                                  % Message length at a time

% Convolution code
% specifies the delay for the encoder's k input bit streams, changing it
% with operation 4 will plot three cases for constraint length 5 7 9, the
% defalut given length is 9.
constlength = 9;                        % constraint length

% Choose any case for individual evalutation or all in case you choose
% operation 4 and run it one by one and each run change the value of it.
if constlength == 5                     
    polynomial = [37 33];               % Generator polynomial for 5
elseif constlength == 7
    polynomial = [171 133];             % Generator polynomial for 7
elseif constlength == 9
    polynomial = [657 435];             % Generator polynomial for 9
end
traceback = 5 * constlength;            
%% Encoding
if operation == 0                       % Uncoded case
    s_enc = s;                          % Enconded data is the same.
elseif operation == 1                   % Repetition code
    s_enc = repelem(s, rep_num);        % Repeating each element by n int.
elseif operation == 2                   % Linear block case
    pol = cyclpoly(n,k);                % Cyclic generator polynomial
    parmat = cyclgen(n,pol);            % Parity-check matrix
    genmat = gen2par(parmat);           % Converted to a generator matrix
    % Encoding the signal with the generator matrix 
    s_enc = encode(s,n,k,'linear/binary',genmat);
else
    % Trellis map generation from the givem polynomial that ecodes the data
    trellis = poly2trellis(constlength, polynomial);
    s_enc = convenc(s, trellis);        % Perfoming the actual encoding
end

n_s_enc = length(s_enc);                % Storing the lenth of encoded data

%% BPSK modulation
% It can be base or pass band, pass band was choosen as not mentiond which
% one to be use, also PRK (phase reversal key) is used; as it's very
% common and guranted to be better facing noise.
fc = 1e6;                               % Carrier frequency
Tb = 1 / fc;                            % Bit period
M = 2;                                  % PSK mod ==> 2: BPSK, 4: QPSK

s_mod = pskmod(s_enc, M, pi);           % NRZ encoder

t_s_mod = (0: Tb: length(s_mod)*Tb-Tb)';% Time of the signal with step Tb
s_tx = s_mod .* cos(2*pi*fc*t_s_mod);   % Carrier multiplication

%% Channel simulation
SNR_dB = 0:2:30;                        % SNR in dB 0 2 4 ... 30
snr_len = length(SNR_dB);               % length of it's array
err_num = zeros(1, snr_len);            % Memory allocation for error nom.

% These will hold the data of probability of error in case of repetition
% codes with n = 3, 5, 11
Pe3 = zeros(1, snr_len);                % M.A. for pe in case of rep = 3
Pe5 = zeros(1, snr_len);                % M.A. for pe in case of rep = 5
Pe11 = zeros(1, snr_len);               % M.A. for pe in case of rep = 11

% These will hold the data of Bit Rate error in case of convlution codes
% with constraint length = 5, 7, 9
BER5 = zeros(1, snr_len);               % M.A. for BER in case of c.l = 5 
BER7 = zeros(1, snr_len);               % M.A. for BER in case of c.l = 7 
BER9 = zeros(1, snr_len);               % M.A. for BER in case of c.l = 9

% This loop is for testing all cases of SNR with different encoding
% techniques
for snr_itr = 1:16
    % Adding AWGN to signal by SNR
    s_tx_noise = awgn(s_tx, SNR_dB(snr_itr));
    
    %% Receiver end
    s_rx = s_tx_noise.*cos(2*pi*fc*t_s_mod);% Carrier multiplication

    s_dem = pskdemod(s_rx, M, pi);          % BPSK demodulation( Vth = 0);

    %% Decoding 
    if operation == 0                       % Without encoding case
        s_dec = s_dem;                      % No decoding is performed
        
    elseif operation == 1                   % Repetition code case
        % This function performs integration(low pass filter), then
        % averages the samples for a certain period, which will be the
        % repetition times. 1 0 1 => 0.67. Weighted average. It also cannot
        % produce 0.5; because we have odd number of repetition all time.
        s_dec = intdump(s_dem, rep_num);    % Integrate and average 
        s_dec(s_dec > 0.5) = 1;             % If sample > .5 ==> 1
        s_dec(s_dec < 0.5) = 0;             % If sample < .5 ==> 0
        
    elseif operation == 2                   % Linear block case
        % Decoding by the decode function, which is simply the reverse
        % operation of encode.
        s_dec = decode(s_dem,n,k,'linear/binary',genmat);
        
        s_dec = s_dec(1: n_sym);            % Removing any null points
    else
        % This function performs the decoding of convlution code by using
        % the trellis map, generted above.
        s_dec = vitdec(s_dem, trellis, traceback, 'trunc', 'hard');
    end
    
    % At each decoding operation test error bits form the original signal
    err_num(snr_itr) = biterr(s_dec, s);
    
    % Calculation probability of error in case of repetition code. with all
    % the given snr; as required 
    if rep_num == 3
        Pe3(snr_itr) = mean(abs(s_dec-s));  % incase of 3
    elseif rep_num == 5
        Pe5(snr_itr) = mean(abs(s_dec-s));  % in case of 5
    elseif rep_num == 11
        Pe11(snr_itr) = mean(abs(s_dec-s)); % In Case of 11
    end
    
end

BER = err_num / n_sym;                      % Assigning BER for general use

if constlength == 5                         % Filling the BER in case of 5
    BER5 = BER;
elseif constlength == 7                     % Filling the BER in case of 7
    BER7 = BER;
else                                        % Filling the BER in case of 9
    BER9 = BER;
end

%% Plottin the comparisions
if operation == 0                               % Uncoded case
    semilogy(SNR_dB, BER, 'y')
    title('BER of non encoding signal');
    axis([0 13 10e-6 10e-1]);
    xlabel('SNR_dB'); 
    ylabel('BER');
    grid on;
    hold on;
  
elseif operation == 1                           % Repetition code
    semilogy(SNR_dB,Pe3, 'r')
    grid on;
    hold on, 
    semilogy(SNR_dB,Pe5, 'g')
    semilogy(SNR_dB,Pe11, 'b')
    title('Pe of repetition codes signal with different n');
    legend('n = 3', 'n = 5', 'n =11');
    axis([0 11 10e-7 10e-2]);
    ylabel('SNR_dB');
    xlabel('Pe');
    
elseif operation == 2                           % Linear block code
    semilogy(SNR_dB, BER, 'c')
    title('BER vs SNR for linear block code');
    axis([0 11 10e-7 10e-1]);
    xlabel('SNR'); 
    ylabel('BER')
    grid on;
    hold on;

elseif operation == 3                           % Convlution code
    semilogy(SNR_dB, BER, 'm')
    title('BER vs SNRn(dB) for linear block code');
    axis([0 11 10e-7 10e-1]);
    legend('NC', 'LPC', 'CC');
    xlabel('SNR'); 
    ylabel('BER')
    grid on;
    hold on;
    
elseif operation == 4
    semilogy(SNR_dB, BER5, 'r')
    hold on;
    semilogy(SNR_dB, BER7, 'g')
    semilogy(SNR_dB, BER9, 'b')
    title('BER for convlution code with different constraint lengths');
    axis([0 5 10e-7 10e-1]);
    legend('constlength=5', 'constlength=7', 'constlength=9');
    xlabel('SNR'); 
    ylabel('BER')
    grid on;
    hold on;
end











