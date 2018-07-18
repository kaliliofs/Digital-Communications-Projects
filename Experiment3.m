close all; clear all; clc;

NumberOfBits = 1e6;             %Digital Signal
SNR_dB = (0:2:30);              %SNR in dB
SNR_watt = 10.^(SNR_dB./10);    %SNR in Watt
%Array Of Binary Bits With Size Of "NumberOfBits"
BinaryBits = randi([0 1],1,NumberOfBits); 
BER11=zeros(1,16);
breakPoint = -1;

%%
%%%%%%%%%%%%%%%% OOK Modulation %%%%%%%%%%%%%%%%%%
signal_OOK = BinaryBits;
breakPoint = -1;

for n = 1 : 16      %The Other 15 SNR values
    
    %Noise = (1/sqrt(snr(n))).*randn(1,NumberOfBits);  %AWGN
    %Rx_sequence = Noise + BinaryBits;  %Receiver
    Rx_sequence = awgn(signal_OOK,SNR_dB(n),'measured');
    %If you use 'measured', then awgn actually measures the signal power first.
    
    Threshold = 0.5.*ones(1,NumberOfBits);
    %If "Rx_sequence" Is Greater Than "Threshold",
     %Then "ReceivedSignal=1" And Vice Versa.
    ReceivedSignal = gt(Rx_sequence,Threshold);
    
    %comparing the transmitted and the received signal
    BitsInError = xor(ReceivedSignal,signal_OOK); 
    NumberOfErrors(n)=sum(BitsInError,2);
    %another method
    [NumberOfErrors11,BER11(n)]=biterr(ReceivedSignal,signal_OOK);
    
    if(NumberOfErrors(n) == 0 && breakPoint == -1)
        breakPoint = SNR_dB(n);
    end  
end
fprintf('The OOK system started to be without error @ SNR = %d\n', breakPoint);

%with the other method
semilogy(SNR_dB,BER11)
xlabel('SNR'); ylabel('BER11')
legend('another way')
grid on;
figure

subplot(3,1,1)
BER_OOK =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER_OOK)
xlabel('SNR'); ylabel('BER_OOK')
grid on;
hold on

subplot(3,1,2)
PE = 0.5.*erfc(0.5.*sqrt(SNR_watt.*2));     
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;

subplot(3,1,3)
BER1 = [BER_OOK' PE'];
semilogy(SNR_dB,BER_OOK,'b', SNR_dB,PE, 'r')
legend('BER_OOK','PE')
xlabel('SNR'); ylabel('BER1')
grid on;



%%
%%%%%%%%%%%%%%%% PRK Modulation %%%%%%%%%%%%%%%%%%
signal_PRK = 2*BinaryBits-1;
breakPoint = -1;

for n = 1 : 16      %The Other 15 SNR values
    
    %Noise = (1/sqrt(snr(n))).*randn(1,NumberOfBits);  %AWGN
    %Rx_sequence = Noise + BinaryBits;  %Receiver
    Rx_sequence = awgn(signal_PRK,SNR_dB(n),'measured');
    %If you use 'measured', then awgn actually measures the signal power first.
    
    Threshold = zeros(1,NumberOfBits);
    %If "Rx_sequence" Is Greater Than "Threshold",
     %Then "ReceivedSignal=1" And Vice Versa.
    ReceivedSignal = gt(Rx_sequence,Threshold);
    
    %comparing the transmitted and the received signal
    BitsInError = gt(2*ReceivedSignal-1,signal_PRK); 
    NumberOfErrors(n)=sum(BitsInError,2);
    
    if(NumberOfErrors(n) == 0 && breakPoint == -1)
        breakPoint = SNR_dB(n);
    end  
end
fprintf('The PRK system started to be without error @ SNR = %d\n', breakPoint);

subplot(3,1,1)
BER_PRK =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER_PRK)
xlabel('SNR'); ylabel('BER_PRK')
grid on;
hold on

subplot(3,1,2)
PE = 0.5*erfc(sqrt(SNR_watt/2));     %probability_of_error
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;

subplot(3,1,3)
BER1 = [BER_PRK' PE'];
semilogy(SNR_dB,BER_PRK,'b', SNR_dB,PE, 'r')
legend('BER_PRK','PE')
xlabel('SNR'); ylabel('BER1')
grid on;


%%
%%%%%%%%%%%%%%%% orthogonal-FSK Modulation %%%%%%%%%%%%%%%%%%
signal_oFSK = zeros(1, 1e6);
for n = 1 : 1e6
    if BinaryBits(n) == 0
        signal_oFSK(n) = 1;
    else
        signal_oFSK(n) = 1i;
    end
end
breakPoint = -1;

for n = 1 : 16      %The Other 15 SNR values
    
    %Noise = (1/sqrt(snr(n))).*randn(1,NumberOfBits);  %AWGN
    %Rx_sequence = Noise + BinaryBits;  %Receiver
    Rx_sequence = awgn(signal_oFSK,SNR_dB(n),'measured');
    %If you use 'measured', then awgn actually measures the signal power first.
    
    Threshold = (1/sqrt(2))+(1/sqrt(2))*i;
    %If "Rx_sequence" Is Greater Than "Threshold",
     %Then "ReceivedSignal=1" And Vice Versa.
    ReceivedSignal = gt(Rx_sequence,Threshold);
    
    %comparing the transmitted and the received signal
    BitsInError = gt(ReceivedSignal,signal_oFSK); 
    NumberOfErrors(n)=sum(BitsInError,2);
    
    if(NumberOfErrors(n) == 0 && breakPoint == -1)
        breakPoint = SNR_dB(n);
    end  
end
fprintf('The o-FSK system started to be without error @ SNR = %d\n', breakPoint);

subplot(3,1,1)
BER_OFSK =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER_OFSK)
xlabel('SNR'); ylabel('BER_OFSK')
grid on;
legend('OOK','PRK','O-FSK')

subplot(3,1,2)
PE = 0.5*erfc(sqrt(SNR_watt));     %probability_of_error
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;

subplot(3,1,3)
BER1 = [BER_OFSK' PE'];
semilogy(SNR_dB,BER_OFSK,'b', SNR_dB,PE, 'r')
legend('BER_OFSK','PE')
xlabel('SNR'); ylabel('BER1')
grid on;



%%          ******* BONUS *******
% OOK-modulation
BER_ASK=zeros(1,16);
h = modem.pammod('M', 2, 'InputType', 'Bit');   %M is M-arr value=2
pdemod = modem.pamdemod(h);         %fn do the ASK
ASK_Tx = modulate(h, BinaryBits);
for n = 1 : 16
    Rx_sequence = awgn(ASK_Tx,SNR_dB(n),'measured');
    ASK_Rx= demodulate(pdemod , Rx_sequence);
    BitsInError = xor(ASK_Rx,BinaryBits); 
    NumberOfErrors(n)=sum(BitsInError,2);
end

BER =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER)
xlabel('SNR'); ylabel('BER_OOK')
grid on;
figure

PE = 0.5.*erfc(0.5.*sqrt(SNR_watt*2));
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;
figure

BER1 = [BER' PE'];
semilogy(SNR_dB,BER,'b', SNR_dB,PE, 'r')
legend('BER','PE')
xlabel('SNR'); ylabel('BER1')
grid on;


%%          ******* BONUS *******
% PRK-modulation
BER_PSK=zeros(1,16);
h = modem.pskmod('M', 2, 'InputType', 'Bit');
pdemod = modem.pskdemod(h);
PSK_Tx = modulate(h, BinaryBits);

for n = 1 : 16
    Rx_sequence = awgn(PSK_Tx,SNR_dB(n),'measured');
    PSK_Rx= demodulate(pdemod , Rx_sequence);
    BitsInError = xor(PSK_Rx,BinaryBits); 
    NumberOfErrors(n)=sum(BitsInError,2);
end
hold on
BER =NumberOfErrors./NumberOfBits ;
semilogy(SNR_dB,BER)
xlabel('SNR'); ylabel('BER_PRK')
grid on;
figure

PE = 0.5*erfc(sqrt(SNR_watt/2));
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;

BER1 = [BER' PE'];
semilogy(SNR_dB,BER,'b', SNR_dB,PE, 'r')
legend('BER','PE')
xlabel('SNR'); ylabel('BER1')
grid on;

%%          ******* BONUS *******
%QAM -- probability of error caclulation
PE = 0.25*(2*erfc(sqrt(2*SNR_watt/5)) + erfc(3*sqrt(2*SNR_watt/5)) + erfc(5*sqrt(2*SNR_watt/5)));
semilogy(SNR_dB,PE)
xlabel('SNR'); ylabel('PE')
grid on;















