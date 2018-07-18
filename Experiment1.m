clear all; clc;

SNR_dB = 0:2:30;                 %SNR in dB
SNR_watt = 10.^(SNR_dB./10);    %SNR in Watt
NumberOfBits = 1e6;              %Digital Signal

%Array Of Binary Bits With Size Of "NumberOfBits"
BinaryBits = randi([0 1], 1,NumberOfBits); 

for n = 1 : 16      %The Other 15 SNR values
    
    %Noise = (1/sqrt(snr(n))).*randn(1,NumberOfBits);  %AWGN
    %Rx_sequence = Noise + BinaryBits;  %Receiver
    Rx_sequence = awgn(BinaryBits,SNR_dB(n),'measured');
    %If you use 'measured', then awgn actually measures the signal power first.
    
    Threshold = 0.5.*ones(1,NumberOfBits);
    %If "Rx_sequence" Is Greater Than "Threshold",
     %Then "ReceivedSignal=1" And Vice Versa.
    ReceivedSignal = gt(Rx_sequence,Threshold);
    
    %comparing the transmitted and the received signal
    BitsInError = xor(ReceivedSignal,BinaryBits); 
   
    NumberOfErrors(n)=sum(BitsInError,2);
    
end

BER =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER)
xlabel('SNR'); ylabel('BER')
grid on;
figure

PE = 0.5.*erfc(0.5.*sqrt(SNR_watt./2));     %probability_of_error
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;
figure

BER1 = [BER' PE'];
semilogy(SNR_dB,BER,'b', SNR_dB,PE, 'r')
legend('BER','PE')
xlabel('SNR'); ylabel('BER1')
grid on;



