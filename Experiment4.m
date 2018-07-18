
%                       Pulse Code Modulation 
%% ====================================================================
clear 
clc


%Quantization error

Amplitude=1;
freq=2;
sampling_freq=4000;     %larger than the nyquest rate -->Fn=4
t=0:(1/sampling_freq):1/freq;
y=Amplitude*cos(2*pi*freq*t);

figure
plot(t,y)
ylabel('Amplitude'); xlabel('time');
legend('Signal before quantization')
figure
%******************
n=[3,4,5,10];
m=2.*n+1;
%*******************
for i=1:length(n)
   y_quantized= double(fi(y,1,m(i),n(i)));
   quantization_error = sum((y_quantized - y).^2);  %qe from the equation
   subplot(2,2,i)
   plot(t,y_quantized,'b')
   ylabel('Amplitude'); xlabel('time');
   legend(['qe= ',num2str(quantization_error)])
   title(['at n= ',num2str(n(i)),'bit'])

end


%% ====================================================================
%signal in binary form

n=input('choose which n you want to encode with ( n[3,4,5,10] ): ');
m=2*n+1;
y_quantized= double(fi(y,1,m,n));
   
   
index=30;
valuesOfQuantization(1:(length(y_quantized)-1)/2)=zeros(1,(length(y_quantized)-1)/2);
valuesOfQuantization(1:(length(y_quantized)-1)/2)=y_quantized(1:(length(y_quantized)-1)/2);
for l=1:3   %shrinking the array to be with single values, not repeted values
    for c=1:length(valuesOfQuantization)
        for j=1+c:length(valuesOfQuantization)
            if(valuesOfQuantization(c)==valuesOfQuantization(j))
                valuesOfQuantization(j:end-1)=valuesOfQuantization(j+1:end);
                j=j-1;
                valuesOfQuantization(end)=index;
                index=index+1;
            end
        end
    end
end
l=0;
c=0;

%only a method used to get the single valued signal
index = find(valuesOfQuantization==30)-1;
valuesOfQuantization2=valuesOfQuantization(1:index);
i=(1:index)';
binary=de2bi(index-i); 

binary_signal=zeros(length(y_quantized),length(binary(1,:)));
for j=1:length(y_quantized)
    for l=1:length(valuesOfQuantization2)
        if valuesOfQuantization2(l)==y_quantized(j)
            binary_signal(j,:)=binary(l,:);
        end
    end
end

j=(1:length(y_quantized))';
disp('    Dec          Binary       ')
disp('   -----   -------------------')
disp([y_quantized', binary_signal])


%% ========================================================================================================
%sampling distorsion

%% ====================================================================
%reconstruction from oversampling

clear 
clc

Fs=1000;
t=0:1/Fs:1;   
y=2*cos(2*pi*5*t);

[B,A] = butter(3,1000/100000,'low');
%it gives the parameters of a lowpass Butterworth filter,
%this lowpass filter give desired cutoff frequency  on 1000 hz from 100000hz, also it's in 3rd order
zero_added_signal=zeros(1,length(y)*10);
for i = 1:length(y)
    zero_added_signal(i*10)=y(i);
end

zero_added_signal(1:9)=[];
%Adding zeros enhances the signal display and don't change
%the spectrum, it changes sampling freq. only
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
plot(t,filtered_signal,'r')
xlabel('time')
ylabel('oversampled signals')

s=fft(filtered_signal);
s=fftshift(s);
fs=Fs;         %***why 100? 
               %must be the same as sampling rate of the original
               %signal above, because we can't change the axis that we
               %plot with, also in upsampling, bandwidth must be decreased
               %not increased
freq=linspace(-fs/2,fs/2,length(s));
figure
plot(freq,abs(s))
xlim([-2 2]);
xlabel('freq')
ylabel('magnitude of oversampled signals')

%% =====================================================================
%construction from minimum sampling
clear 
clc

figure
Fs=10;      % minimum sampling Fn = 2*Fm
t=0:1/Fs:1;    
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.1,'low');
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);   %upsampling -> decreases BW
end

zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
plot(t,filtered_signal,'r')
xlabel('time')
ylabel('minimum sampling signals')

s=fft(filtered_signal);
s=fftshift(s);
fs=Fs;         %***why 100? 
               %must be the same as sampling rate of the original
               %signal above, because we can't change the axis that we
               %plot with, also in upsampling, bandwidth must be decreased
               %not increased
freq=linspace(-fs/2,fs/2,length(s));
figure
plot(freq,abs(s))
xlim([-2 2]);
xlabel('freq')
ylabel('magnitude of minimum sampled signals')

%% ====================================================================
%construction from undersampling sampling
clear 
clc

figure
Fs=10/2;        % undersampling sampling Fn < 2*Fm
t=0:1/Fs:1;    
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.2,'low');

zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);   %upsampling -> decreases BW
end

zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
plot(t,filtered_signal,'r')
xlabel('time')
ylabel('undersampling sampling signals')

s=fft(filtered_signal);
s=fftshift(s);
fs=Fs;         %***why 50?
               %must be the same as sampling rate of the original
               %signal above, because we can't change the axis that we
               %plot with, also in upsampling, bandwidth must be decreased
               %not increased
freq=linspace(-fs/2,fs/2,length(s));
figure
plot(freq,abs(s))
xlim([-2 2]);
xlabel('freq')
ylabel('magnitude of undersampled signals')


%% ========================================================================================
%-----------------------------------------------------------------------------------------

%bouns2
clear 
clc

%Quantization with "quantize"

Amplitude=1;
freq=2;
sampling_freq=4000;     %larger than the nyquest rate -->Fn=4
t=0:(1/sampling_freq):2*1/freq;
y=Amplitude*cos(2*pi*freq*t);

%******************
n=[3,4,5,10];
m=2.*n+1;
%*******************
for i=1:length(n)
   q = quantizer('fixed', [m(i) n(i)]); 
   %a quantizer object with properties set to its inputs.
   
   y_quantized=quantize(q, y);
   quantization_error = sum((y_quantized - y).^2);  %--> quantization error from eqn
   subplot(2,2,i)
   plot(t,y_quantized,'b')
   legend(['qe= ',num2str(quantization_error)])
   ylabel('Amplitude'); xlabel('time');
   title(['at n= ',num2str(n(i)),'bit'])

end

                % ====================================

%signal in binary form

n=input('choose which n you want to encode with ( n[3,4,5,10] ): ');
m=2*n+1;
y_quantized= double(fi(y,1,m,n));
   
   
index=30;
valuesOfQuantization(1:(length(y_quantized)-1)/2)=zeros(1,(length(y_quantized)-1)/2);
valuesOfQuantization(1:(length(y_quantized)-1)/2)=y_quantized(1:(length(y_quantized)-1)/2);
for l=1:4   %shrinking the array to be with single values, not repeted values
    for c=1:length(valuesOfQuantization)
        for j=1+c:length(valuesOfQuantization)
            if(valuesOfQuantization(c)==valuesOfQuantization(j))
                valuesOfQuantization(j:end-1)=valuesOfQuantization(j+1:end);
                j=j-1;
                valuesOfQuantization(end)=index;
                index=index+1;
            end
        end
    end
end
l=0;
c=0;

%only a method used to get the single valued signal
index = find(valuesOfQuantization==30)-1;
valuesOfQuantization2=valuesOfQuantization(1:index);
i=(1:index)';
binary=de2bi(index-i); 

binary_signal=zeros(length(y_quantized),length(binary(1,:)));
for j=1:length(y_quantized)
    for l=1:length(valuesOfQuantization2)
        if valuesOfQuantization2(l)==y_quantized(j)
            binary_signal(j,:)=binary(l,:);
        end
    end
end

j=(1:length(y_quantized))';
disp('    Dec          Binary       ')
disp('   -----   -------------------')
disp([y_quantized', binary_signal])

%% --------------------------------------------------------------------------------------
%bouns3

clear
clc
%Non uniform Quantization with "compand"

Amplitude=1;
freq=2;
sampling_freq=4000;     %larger than the nyquest rate -->Fn=4
t=0:(1/sampling_freq):2*1/freq;
y=1*cos(2*pi*freq*t);

%******************
n=[3,4,5,10];
m=2.*n+1;
compressed = compand(y,255,max(y),'mu/compressor');
%-Non-uniform quantization is achieved by, first passing the input signal through a “compressor”. 
  %The output of the compressor is then passed through a uniform quantizer.
%-The combined effect of the compressor and the uniform quantizer is that of a non-uniform quantizer.
%At the receiver the signal is restored to its original form by using an expander. 

%*******************

for i=1:length(n)
   q = quantizer('fixed', [m(i) n(i)]);
   
   y_quantized=quantize(q, compressed);     %Non uniform quantization
   quantization_error = sum((y_quantized - y).^2);  %--> quantization error from eqn
   subplot(2,2,i)
   plot(t,y_quantized,'b')
   legend(['qe= ',num2str(quantization_error)])
   ylabel('Amplitude'); xlabel('time');
   title(['at n= ',num2str(n(i)),'bit'])


end


                % ====================================

%signal in binary form

n=input('choose which n you want to encode with ( n[3,4,5,10] ): ');
m=2*n+1;
y_quantized= double(fi(y,1,m,n));
   
   
index=30;
valuesOfQuantization(1:(length(y_quantized)-1)/2)=zeros(1,(length(y_quantized)-1)/2);
valuesOfQuantization(1:(length(y_quantized)-1)/2)=y_quantized(1:(length(y_quantized)-1)/2);
for l=1:4   %shrinking the array to be with single values, not repeted values
    for c=1:length(valuesOfQuantization)
        for j=1+c:length(valuesOfQuantization)
            if(valuesOfQuantization(c)==valuesOfQuantization(j))
                valuesOfQuantization(j:end-1)=valuesOfQuantization(j+1:end);
                j=j-1;
                valuesOfQuantization(end)=index;
                index=index+1;
            end
        end
    end
end
l=0;
c=0;

%only a method used to get the single valued signal
index = find(valuesOfQuantization==30)-1;
valuesOfQuantization2=valuesOfQuantization(1:index);
i=(1:index)';
binary=de2bi(index-i); 

binary_signal=zeros(length(y_quantized),length(binary(1,:)));
for j=1:length(y_quantized)
    for l=1:length(valuesOfQuantization2)
        if valuesOfQuantization2(l)==y_quantized(j)
            binary_signal(j,:)=binary(l,:);
        end
    end
end

j=(1:length(y_quantized))';
disp('    Dec          Binary       ')
disp('   -----   -------------------')
disp([y_quantized', binary_signal])




