clc;
clear all;
close all;
%%
%Μέρος 1ο - Ψυχοακουστικό Μοντέλο

%Βήμα 1.0

[temp,Fs] =audioread('music-dsp19.wav');
%conversion to monophonic
music = 0.5*(temp(:,1)+temp(:,2));
%signal normalization
music = music./max(abs(music(:)));

t = (0:length(music)-1)/Fs;
figure(1);
plot(t,music);
title('Initial Music Signal');
xlabel('Time (sec)');


N = 512;
wHann = hanning(N);

% Framing and windowing the music signal
sigFramed = buffer(music, N,0);
x = (diag(sparse(wHann)))*sigFramed;

figure(2);
t2 = (512*634+1:512*635)/Fs;
plot(t2,x(:,635));
title('Middle (635th) window of initial signal');
xlabel('time(sec)');


%Βήμα 1.1

p = 90.302 + 10*log10(abs(fft(x)).^2);
p = p(1:N/2,:);

%Βήμα 1.2

%Tone Maskers
Ptm = zeros(N/2,size(p,2));
for i=1:size(p,2)
    Ptm(:,i) = findToneMaskers(p(:,i));
end

%Noise Maskers
f = (1:N/2)*(Fs/N);
b = hz2bark(f);
Pnm = zeros(N/2,size(p,2));
for i = 1:size(p,2)
    Pnm(:,i) = findNoiseMaskers(p(:,i),Ptm(:,i),b);
end


figure(3);
plot(f,p(:,635));
hold on;
stem(f,Ptm(:,635),'-x');
stem(f,Pnm(:,635),'-x');
title('Spectral Density & Masks of 635th window');
xlabel('Frequency(Hz)');
legend({'Spectral Density','Tone Maskers','Noise Maskers'});
hold off;



figure(4);
plot(b,p(:,635));
hold on;
stem(b,Ptm(:,635),'-x');
stem(b,Pnm(:,635),'-x');
title('Spectral Density & Masks of 635th window');
xlabel('Frequency(Bark)');
legend({'Spectral Density','Tone Maskers','Noise Maskers'});
hold off;


%Βήμα 1.3

Tq = 3.64*(f./1000).^(-0.8)-6.5*exp(-0.6*(f./1000-3.3).^2)+10.^(-3)*(f./1000).^4;
for i=1:size(p,2)
    [Ptm(:,i),Pnm(:,i)] = checkMaskers(Ptm(:,i).',Pnm(:,i).',Tq,b);
end

figure(5);
plot(f,p(:,635));
hold on;
stem(f,Ptm(:,635),'-x');
stem(f,Pnm(:,635),'-x');
title('Spectral Density & Final Masks of 635th window');
xlabel('Frequency(Hz)');
legend({'Spectral Density','Tone Maskers','Noise Maskers'});
hold off;


figure(6);
plot(b,p(:,635));
hold on;
stem(b,Ptm(:,635),'-x');
stem(b,Pnm(:,635),'-x');
title('Spectral Density & Final Masks of 635th window');
xlabel('Frequency(Bark)');
legend({'Spectral Density','Tone Maskers','Noise Maskers'});
hold off;


%Βήμα 1.4

Ttm = zeros(N/2,N/2,size(p,2));
Tnm = zeros(N/2,N/2,size(p,2));
for i=1:size(p,2)
    [Ttm(:,:,i),Tnm(:,:,i)] = IndividualMaskingThresholds(Ptm(:,i),Pnm(:,i),b);
end


%Βήμα 1.5

Tg = zeros(N/2,size(p,2));
for i=1:size(p,2)
    Tg(:,i) = GlobalMaskingThreshold(Ttm(:,:,i),Tnm(:,:,i),Tq.');
end

figure(7);
plot(f,p(:,635));
hold on;
stem(f,Ptm(:,635),'-x');
stem(f,Pnm(:,635),'-x');
plot(f,Tg(:,635));
title('Global Masking Threshold of 635th window');
xlabel('Frequency(Hz)');
legend({'Spectral Density','Tone Maskers','Noise Maskers','Global Masking Threshold'});
hold off;



figure(8);
plot(b,p(:,635));
hold on;
stem(b,Ptm(:,635),'-x');
stem(b,Pnm(:,635),'-x');
plot(b,Tg(:,635));
title('Global Masking Threshold of 635th window');
xlabel('Frequency(Bark)');
legend({'Spectral Density','Tone Maskers','Noise Maskers','Global Masking Threshold'});
hold off;

%%
%Μέρος 2ο - Χρονο-Συχνοτική Ανάλυση με Συστοιχία Ζωνοπερατών Φίλτρων

%Βήμα 2.0

[h,g]=Filterbank(32);

%Βήμα 2.1-2.3
%Adaptive Quantization

for i=1:size(sigFramed,2)
     [synthesized1(:,i),bits1(i)] =QuantizationSinthesis(sigFramed(:,i),h,g,Tg(:,i),32,1);
end

totalbits1 = sum(bits1);

%convolution with filters adds 64 = 2M samples delay at the start and the end of each widnow
synthesized1 = synthesized1(64:575,:);
k=1;
new1 = zeros(1,size(synthesized1,1)*size(synthesized1,2));
for j=1:size(synthesized1,2)
    for i=1:size(synthesized1,1)
    new1(k) = synthesized1(i,j);
    k=k+1;
    end
end
new1 = new1(1:length(music));
audiowrite('C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ-Εργαστήριο 2019\Lab2\music_adaptive.wav',new1,44100);

t = (0:length(new1)-1)/Fs;
figure(9);
plot(t,new1);
title('Compressed Music Signal using Psychoacoustic Model');
xlabel('Time (sec)');
ylim([-1 1]);


%Non-Adaptive Quantization

for i=1:size(sigFramed,2)
     [synthesized2(:,i),bits2(i)] =QuantizationSinthesis(sigFramed(:,i),h,g,Tg(:,i),32,0);
end

totalbits2 = sum(bits2);

%convolution with filters adds 64 = 2M samples delay at the start and the end of each widnow
synthesized2 = synthesized2(64:575,:);
k=1;
new2 = zeros(1,size(synthesized2,1)*size(synthesized2,2));
for j=1:size(synthesized2,2)
    for i=1:size(synthesized2,1)
    new2(k) = synthesized2(i,j);
    k=k+1;
    end
end
new2 = new2(1:length(music));
audiowrite('C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ-Εργαστήριο 2019\Lab2\music_8bit.wav',new2,44100);

t = (0:length(new2)-1)/Fs;
figure(10);
plot(t,new2);
title('Compressed Music Signal using 8 bits per sample');
xlabel('Time (sec)');
ylim([-1 1]);


%comparison between the 2 methods

%compression ratio for adaptive quantization
compression1 = totalbits1 /(16*length(music));

%compression ratio for non-adaptive quantization
compression2 = totalbits2 /(16*length(music));

%mse for adaptive quantization
mse1 = sum((music-new1(1,1:650415).').^2)/length(music);

%mse for non-adaptive quantization
mse2 = sum((music-new2(1,1:650415).').^2)/length(music);

t = (0:length(new1)-1)/Fs;
figure(11);
plot(t,(music-new1(1,1:650415).').^2);
title('Square Error for adaptive Quantization');
xlabel('Time (sec)');


t = (0:length(new2)-1)/Fs;
figure(12);
plot(t,(music-new2(1,1:650415).').^2);
title('Square Error for non-adaptive Quantization');
xlabel('Time (sec)');

