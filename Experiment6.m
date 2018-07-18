clear all;
close all;
%Number of bits to be modulated
NumberOfBits=10;

%Generate some bits randomly
RandomNumber=rand(1,NumberOfBits)>0.5;

%Bit rate in bits/second
BitRate=4;

%Assuming waveform is represented by 10 samples so sampling
%frequency is 10 times the bit rate
SamplingFrequncy=10*BitRate;

%Define arrays for output waves
NRZ_out=[];
NRZI_out=[];
RZ_out=[];
AMI_out=[];
Manchester_out=[];
MLT3_out=[];
 
%Peak voltage (+vp)
PeakVoltage=5;

%NRZ Modulation 
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage];
 elseif RandomNumber(index)==0
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-PeakVoltage)];
 end
end

%NRZ-Inverted Modulation
%OneFlag is a flag used to indicate the last "One" state (positive/negative)
OneFlag = 1; %Initial value from +vp
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1
 OneFlag = -1*OneFlag; %Invert the "One" state
 NRZI_out=[NRZI_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage*OneFlag];
 elseif RandomNumber(index)==0
 NRZI_out=[NRZI_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage*OneFlag];
 end
end

%RZ Modulation
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*PeakVoltage];
 elseif RandomNumber(index)==0
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*(-PeakVoltage)];
 end
end

%AMI Modulation
%OneFlag is a flag used to indicate the last "One" state (positive/negative)
OneFlag = 1; %Initial value from +vp
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1
 OneFlag = -1*OneFlag; %Invert the "One" state
 AMI_out=[AMI_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage*OneFlag];
 elseif RandomNumber(index)==0
 AMI_out=[AMI_out [0 0 0 0 0 0 0 0 0 0]];
 end
end


%Manchester Modulation
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1
 Manchester_out=[Manchester_out [1 1 1 1 1 -1 -1 -1 -1 -1]*PeakVoltage];
 elseif RandomNumber(index)==0
 Manchester_out=[Manchester_out [1 1 1 1 1 -1 -1 -1 -1 -1]*(-PeakVoltage)];
 end
end

%MLT-3 Modulation
%Level is a vector to indicate the last output level (-1,0,1)
Level = [0 1 0 -1]; 
i = 1;
for index=1:size(RandomNumber,2)
 if RandomNumber(index)==1    
 if (i < 4)
    i = i+1;
 else
    i = 1;
 end
 MLT3_out=[MLT3_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage*Level(i)];
 elseif RandomNumber(index)==0
 MLT3_out=[MLT3_out [1 1 1 1 1 1 1 1 1 1]*PeakVoltage*Level(i)];
 end
end

%Darwing Line code modulated waves
%Draw NRZ
subplot(3,2,1)
stem((1:NumberOfBits*10)/10,NRZ_out);
xlabel('Bits');
ylabel('Voltage');
title('NRZ Modulated wave');

%Draw NRZI
subplot(3,2,2)
stem((1:NumberOfBits*10)/10,NRZI_out);
xlabel('Bits');
ylabel('Voltage');
title('NRZ-Inverted Modulated wave');

%Draw RZ
subplot(3,2,3)
stem((1:NumberOfBits*10)/10,RZ_out);
xlabel('Bits');
ylabel('Voltage');
title('RZ Modulated wave');

%Draw AMI
subplot(3,2,4)
stem((1:NumberOfBits*10)/10,AMI_out);
xlabel('Bits');
ylabel('Voltage');
title('AMI Modulated wave');

%Draw Manchester
subplot(3,2,5)
stem((1:NumberOfBits*10)/10,Manchester_out);
xlabel('Bits');
ylabel('Voltage');
title('Manchester Modulated wave');

%Draw MLT-3
subplot(3,2,6)
stem((1:NumberOfBits*10)/10,MLT3_out);
xlabel('Bits');
ylabel('Voltage');
title('MLT-3 Modulated wave');

%Darwing Welch PSD of the modulated signals
%Draw NRZ
h = spectrum.welch;
Hpsd=psd(h,NRZ_out,'Fs',SamplingFrequncy);
figure;
subplot(3,2,1)
handle=plot(Hpsd);
set(handle,'LineWidth',2.5,'Color','r')

%Draw NRZ-Inverted
h = spectrum.welch;
Hpsd=psd(h,NRZI_out,'Fs',SamplingFrequncy);
subplot(3,2,2)
handle=plot(Hpsd);
set(handle,'LineWidth',2.5,'Color','r')

%Draw RZ
h = spectrum.welch;
Hpsd=psd(h,RZ_out,'Fs',SamplingFrequncy);
subplot(3,2,3)
handle=plot(Hpsd);
set(handle,'LineWidth',2.5,'Color','b')

%Draw AMI
h = spectrum.welch;
Hpsd=psd(h,AMI_out,'Fs',SamplingFrequncy);
subplot(3,2,4)
handle=plot(Hpsd)
set(handle,'LineWidth',2.5,'Color','b')

%Draw Manchester
h = spectrum.welch;
Hpsd=psd(h,Manchester_out,'Fs',SamplingFrequncy);
subplot(3,2,5)
handle=plot(Hpsd)
set(handle,'LineWidth',2.5,'Color','k')

%Draw MLT-3
h = spectrum.welch;
Hpsd=psd(h,MLT3_out,'Fs',SamplingFrequncy);
subplot(3,2,6)
handle=plot(Hpsd)
set(handle,'LineWidth',2.5,'Color','k')