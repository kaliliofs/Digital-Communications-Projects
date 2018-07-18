%% ========================================================================
% Declaration & initialization 
clear all;  clc;
SNRRange=0:2:30;
m=10;                             %number of samples that represents waveform
NumberOfBits=1e6;                 
filter = ones(10,1);
breakPoint = -1;                  %error becomes zero @ SNR = breakPoint
ber = zeros(1,16);                % from 1 to 16 (range of SNR)
data = randi([0 1],NumberOfBits,1); %Random generation of binary bits 

%% ========================================================================
% Testing the system
for n=1:16
    SNR = SNRRange(n);     
    fprintf('@SNR = %d\n',SNR);
    dataTransmitted = reshape((data*ones(1,m))',NumberOfBits*m,1);
    dataRecieved    = awgn(dataTransmitted,SNR,'measured');
    %MF filter
    MFdataFiltered  = conv(dataRecieved,filter);
    MFdataSampled   = MFdataFiltered(m:10:end);
    MFDifference    = xor((MFdataSampled >= 5),data); %xor between Received signal after comparator with bits
    Num_of_Err      = sum(MFDifference(:));     % (:) to make it column vector
    ber(n)          = Num_of_Err/NumberOfBits;
    %correlator
    corrData        = reshape(dataRecieved,10,NumberOfBits)' * filter;
    CorrDifference  = xor((corrData >= 5),data);    %xor between Received signal after comparator with bits
    if(sum(MFDifference  ~= CorrDifference) == 0)
        fprintf('correlator & MF outputs are the same\n');
    end
    %breakPoint detection
    if(Num_of_Err == 0 && breakPoint == -1)
        breakPoint = SNR;
    end    
end
%% ========================================================================
% plotting ber Vs SNR
tansmisttedPower = (1/NumberOfBits) * sum(data.^2);
display(tansmisttedPower);
fprintf('the system started to be without error @ SNR = %d\n', breakPoint);
semilogy(SNRRange, ber);
xlabel('SNR'); 
ylabel('BER');
grid on;