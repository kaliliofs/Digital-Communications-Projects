

%%
%         ******** sin wave **********
close all;
clear all;  clc;

Fm=500;
Fs=200000;     %100 samples
z=0:(1/Fs):1;
x=3*sin(2*pi*Fm*z);
%===========================
TsScale=20;
k=0:(1/(Fs/TsScale)):1;
delta_signal=zeros(1,length(k));
n=0;
StepSize= 0.3;

%*********************
t=1;
delta_signal(t)=n;
%*********************

for t=2:(Fs/TsScale);
    if (x((t*TsScale))>=n)
        n=n+StepSize;
    elseif(x((t*TsScale))<n)
        n=n-StepSize;
    end
    
    delta_signal(t)=n;
end
B=[0];
figure;
subplot(3,1,1);a=stairs(k,delta_signal)      %Plot Modulated signal
xlabel('time')
ylabel('amplitude signal')
%plot(z,delta_signal);
hold on
plot(z,x); %Plot Original Signal on the same plot of Modulated signal
grid
xlim([0 .005]); %limit X axis
ylim([-3.5 3.5]); %limit Y axis
%Get the Output modulated signal 
for t=1:(Fs/TsScale);
    if (x((t*TsScale))>delta_signal(t))
       B(t+1)=1;
    elseif (x((t*TsScale))==delta_signal(t))
        B(t+1)=1;
    else
        B(t+1)=0;
    end 
end
 subplot(3,1,2);stairs(k,B,'linewidth',3);
 xlabel('time')
 ylabel('delta modulation output')
  xlim([0 .005]);
ylim([0 1]);
grid

%       *********** Filtering The Delta Signal*********
[B,A] = butter(3,10000/100000,'low');
filtered_signal = filter(B,A,delta_signal);
subplot(3,1,3)
plot(k,filtered_signal,'r')
hold on
plot(z,x,'b');
legend('After Filter','original signal')
xlim([0 0.005]);
ylim([-3.5 3.5]);

%        **********Square Error********
for t=1:(Fs/TsScale);
    sqrtErr=power(x((t*TsScale))-filtered_signal(t),2);
end

fprintf('Square Error = %d\n',sqrtErr)


%%
%         ******** DC voltage **********
close all;
clear all;  clc;

Fs=200000;     %100 samples
z=0:(1/Fs):1;
x=ones(1,Fs+1);
%==========================
TsScale=20;
k=0:(1/(Fs/TsScale)):1;
delta_signal=zeros(1,length(k));
n=0;
StepSize= 0.3;

%*********************
t=1;
delta_signal(t)=n;
%*********************

for t=2:(Fs/TsScale);
    if (x((t*TsScale))>=n)
        n=n+StepSize;
    elseif(x((t*TsScale))<n)
        n=n-StepSize;
    end
    
    delta_signal(t)=n;
end
B=[0];
figure;
subplot(3,1,1);a=stairs(k,delta_signal)      %Plot Modulated signal
xlabel('Time')
ylabel('Amp')
%plot(z,delta_signal);
hold on
plot(z,x); %Plot Original Signal on the same plot of Modulated signal
grid
xlim([0 .005]); %limit X axis
ylim([-3.5 3.5]); %limit Y axis
%Get the Output modulated signal 
for t=1:(Fs/TsScale);
    if (x((t*TsScale))>delta_signal(t))
       B(t+1)=1;
    elseif (x((t*TsScale))==delta_signal(t))
        B(t+1)=1;
    else
        B(t+1)=0;
    end 
end
 subplot(3,1,2);stairs(k,B,'linewidth',3);
 xlabel('Time')
 ylabel('Bits')
  xlim([0 .005]);
ylim([0 1]);
grid


%       *********** Filtering The Delta Signal*********
[B,A] = butter(3,10000/100000,'low');
filtered_signal = filter(B,A,delta_signal);
subplot(3,1,3)
plot(k,filtered_signal,'r')
hold on
plot(z,x,'b');
legend('After Filter','original signal')
xlim([0 0.005]);
ylim([-3.5 3.5]);

%        **********Square Error********
for t=1:(Fs/TsScale);
    sqrtErr=power(x((t*TsScale))-filtered_signal(t),2);
end

fprintf('Square Error = %d\n',sqrtErr)

%%
%         ******** square wave **********
close all;
clear all;  clc;

Fs=500;
z=0:(1/Fs):1;
x = ones(1,Fs+1);
x(0.25*Fs:0.75*Fs)=0;
%====================
TsScale=20;
k=0:(1/(Fs/TsScale)):1;
delta_signal=zeros(1,length(k));
n=0;
StepSize= 0.4;

%*********************
t=1;
delta_signal(t)=n;
%*********************

for t=2:(Fs/TsScale);
    if (x((t*TsScale))>n)
        n=n+StepSize;
    elseif(x((t*TsScale))<n)
        n=n-StepSize;
    end
    
    delta_signal(t)=n;
end
figure
subplot(3,1,1);
stairs(k,delta_signal)
hold on
plot(z,x);
grid
xlim([0 1]);
ylim([-3.5 3.5]);
B=[0];
for t=1:(Fs/TsScale);
    if (x((t*TsScale))>delta_signal(t))
       B(t+1)=1;
    elseif (x((t*TsScale))==delta_signal(t))
        B(t+1)=1;
    else
        B(t+1)=0;
    end 
end

subplot(3,1,2);
stairs(k,B,'linewidth',3);
xlabel('time')
ylabel('delta modulation bits')
xlim([0 1]);
ylim([0 1]);
grid
%       *********** Filtering The Delta Signal*********
[B,A] = butter(1,10000/100000,'low');
filtered_signal = 2*filter(B,A,delta_signal);
subplot(3,1,3);
plot(k,filtered_signal,'r')
hold on
plot(z,x,'b');
legend('After Filter','original signal')
xlim([0 1]);
ylim([-3.5 3.5]);

%        **********Square Error********
for t=1:(Fs/TsScale);
    sqrtErr=power(x((t*TsScale))-filtered_signal(t),2);
end

fprintf('Square Error = %d\n',sqrtErr)

