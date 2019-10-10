clc;
clear all;

%%
%Μέρος 1ο
%Ερώτημα 1.1
n=0:999;
d0= sin(0.7217.*n)+sin(1.0247.*n);
d1= sin(0.5346.*n)+sin(0.9273.*n);
d2= sin(0.5346.*n)+sin(1.0247.*n);
d3= sin(0.5346.*n)+sin(1.1328.*n);
d4= sin(0.5906.*n)+sin(0.9273.*n);
d5= sin(0.5906.*n)+sin(1.0247.*n);
d6= sin(0.5906.*n)+sin(1.1328.*n);
d7= sin(0.6535.*n)+sin(0.9273.*n);
d8= sin(0.6535.*n)+sin(1.0247.*n);
d9= sin(0.6535.*n)+sin(1.1328.*n);

%Ερώτημα 1.2
f4=fftshift(fft(d4));
f6=fftshift(fft(d6));

m4=abs(f4);
m6=abs(f6);


f = (-length(m4)/2:length(m4)/2-1)/length(m4);
figure(1);
subplot(2,1,1);
plot(f,m4);
title('|D4[k]|');
xlabel('Frequency (Hz)');
ylabel('|D4[k]|');

f = (-length(m6)/2:length(m6)/2-1)/length(m6);
subplot(2,1,2);
plot(f,m6);
title('|D6[k]|');
xlabel('Frequency (Hz)');
ylabel('|D6[k]|');

%Ερώτημα 1.3
%’θροισμα ψηφίων ΑΜ=03115035+03115059=06230094

z = zeros(1,100);
signal=[d0 z d6 z d2 z d3 z d0 z d0 z d9 z d4];

wavwrite(signal,8192,'C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ- εργαστηριο 1-2019\dsp19_lab1_Data\tone_sequence.wav');


%Ερώτημα 1.4
%i)

sigFramed = buffer(signal,1000,-100);
rectft = fft(sigFramed);
figure(2);
subplot(2,1,1);
t=(0:length(rectft)-1)/length(rectft);
absrectft = abs(rectft);
plot(t,absrectft);
title('Rectangular window');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
%ii)
% A hamming window is of N=1000 length is chosen, skipping the 100 zeros
% after each tone of the signal
winLen = 1000;
winOverlap = -100;
wHamm = hamming(winLen);

% Framing and windowing the signal
sigFramed = buffer(signal, winLen, winOverlap);
sigWindowed = (diag(sparse(wHamm))) * sigFramed;
hamdft = fft(sigWindowed);
abshamdft = abs(hamdft);
t=(0:length(hamdft)-1)/length(hamdft);
subplot(2,1,2);
plot(t,abshamdft);
title('Hamming window');
xlabel('Frequency(Hz)');
ylabel('Amplitude');

%Ερώτημα 1.5

A = sum(abshamdft(:,:),2);
[pks,locs] = findpeaks(A);

k=1;
for i=1:ceil((size(pks,1))/2)
    if pks(i)>100
        B(k) = locs(i);
        k = k+1;
    end
end
fb = 2*pi.*B./1000;
res = [B; fb];
disp(res);
%Ερώτημα 1.6
ttdecode(signal);

f = fftshift(fft(signal));
e = f.*conj(f);
f = (-length(e)/2:length(e)/2-1)/length(e);

figure(3);
subplot(1,1,1);
plot(f,e);
title('E[k]');
xlabel('Frequency (Hz)');
ylabel('E[k]');

%Ερώτημα 1.7
load('C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ- εργαστηριο 1-2019\dsp19_lab1_Data\my_touchtones.mat','hardSig','easySig');
ttdecode(easySig);
ttdecode(hardSig);

%%
%Μέρος 2ο
%Μέρος 2.1

Fs=1000;

%α)
t = 0:1/Fs:2-1/Fs;
v = randn([1 length(t)]);
x = 2*cos(t.*2*pi*70) + 3*sin(t.*2*pi*140) +v.*0.15;
figure(4);
plot(t,x);
title('Discrete Signal x[n] with Sampling Frequency 1000 Hz');
xlabel('Time (sec)');

%β)
[s,f,t] = spectrogram(x,40,20,[],1000);
s = abs(s);
figure(5);
surf(t,f,s);
title('Short Time Fourier Transform Of Signal : x[n] using 40msec Window and 20 msec overlap');
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
 
%γ)
nt=0.001:0.001:2;
[scale,Fd]=wavescales('morl',Fs);   

DT_CWT=cwtft({x,0.001},'scales',scale,'wavelet','morl'); 
dsig=DT_CWT.cfs;    

figure(6);
surf(nt,Fd,abs(dsig),'edgecolor','none');
title('DT-CWT Of Signal x[n]');
xlabel('Time')
ylabel('Frequency')
zlabel('Signal amplitude')

%Μέρος 2.2
%α)
t = 0:1/Fs:2-1/Fs;
v = randn([1 length(t)]);
x = 1.7*cos(2*pi*90.*t) + 0.15.*v;
x(625)= x(625)+1.7;
x(800)= x(800)+1.7;
figure(7);
plot(t,x);
title('Discrete Signal x[n] with Sampling Frequency 1000 Hz');
xlabel('Time (sec)');

%β)
[s,f,t] = spectrogram(x,40,20,[],1000);
s = abs(s);
figure(8);
contour(t,f,s);
title('Short Time Fourier Transform Of Signal : x[n] using 40msec Window and 20 msec overlap ');
xlabel('Time (sec)');
ylabel('Frequency (Hz)');

%γ)
[scale,Fd]=wavescales('morl',Fs);   
DT_CWT=cwtft({x,0.001},'scales',scale,'wavelet','morl');    
dsig=DT_CWT.cfs;

figure(9);
contour(nt,Fd,abs(dsig));
title('DT-CWT of x[n]');
xlabel('Time(sec)');
ylabel('Frequency (Hz)');

%%

%Μέρος 3ο

%Ερώτημα 3.1
% A sentence stored in SPHERE format is read. The sentence is :
% 'Αλλά και για να επέμβει χρειάζεται συγκεκριμένη καταγγελία'. A function of voicebox is used.
speechSignal =audioread('C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ- εργαστηριο 1-2019\dsp19_lab1_Data\speech_utterance.wav');
Fs=16000;
% A hamming window is chosen
winLen = 0.02*Fs;
winOverlap = 0.02*Fs-1;
wHamm = hamming(winLen);

% Framing and windowing the signal without for loops.
sigFramed = buffer(speechSignal, winLen, winOverlap, 'nodelay');
sigWindowed = (diag(sparse(wHamm))) * sigFramed;
sigWindowed2=sigWindowed.^2;
% Short-Time Energy calculation
energyST = sum(sigWindowed2,1);
% Time in seconds, for the graphs
% The signals in the graphs are normalized
t = (0:length(speechSignal)-1)/Fs;
figure(10);
subplot(1,1,1);
plot(t, speechSignal/max(abs(speechSignal)));
title('speech: Αλλά και για να επέμβει χρειάζεται συγκεκριμένη καταγγελία ');
xlims = get(gca,'Xlim');
hold on;

% Short-Time energy is delayed due to lowpass filtering. This delay is
% compensated for the graph.
delay = (winLen - 1)/2;
plot(t(delay+1:end - delay), energyST/max(energyST), 'r');
xlim(xlims);
xlabel('Time (sec)');


sigdif = sign(speechSignal(2:end))-sign(speechSignal(1:end-1));

% No change of sign at the beginning
sigdif = [0; sigdif];
sigdifFramed = buffer(sigdif, winLen, winOverlap,'nodelay');
sigdifWindowed = diag(sparse(wHamm)) * abs(sigdifFramed);
zcr = sum(sigdifWindowed, 1)/(2*winLen);

t = (0:length(speechSignal)-1)/Fs;
delay = (winLen - 1)/2;
plot(t(delay+1:end-delay), zcr/max(zcr), 'g');
legend({'Speech Utterance','Short-Time Energy','Zero Crossing Rate'});

hold off;


%Ερώτημα 3.2
[speechSignal2,Fs] =audioread('C:\Users\Ioanna\Documents\ΣΗΜΜΥ\Σ\ΨΕΣ\ΨΕΣ- εργαστηριο 1-2019\dsp19_lab1_Data\music_cut.wav');
winLen = 0.02*Fs;
winOverlap = 0.02*Fs-1;
wHamm = hamming(winLen);
sigFramed2 = buffer(speechSignal2, winLen, winOverlap, 'nodelay');
sigWindowed3 = diag(sparse(wHamm)) * sigFramed2;
sigWindowed4=sigWindowed3.^2;
% Short-Time Energy calculation
energyST1 = sum(sigWindowed4,1);
% Time in seconds, for the graphs
t = (0:length(speechSignal2)-1)/Fs;



figure(11);
subplot(1,1,1);
plot(t, speechSignal2/max(abs(speechSignal2)));
title('speech: Music cut');
xlabel('Time (sec)');
xlims = get(gca,'Xlim');
hold on;

% Short-Time energy is delayed due to lowpass filtering. This delay is
% compensated for the graph.
delay = (winLen - 1)/2;
plot(t(delay+1:end - delay), energyST1/max(energyST1), 'r');
xlim(xlims);
xlabel('Time (sec)');



sigdif2 = sign(speechSignal2(2:end))-sign(speechSignal2(1:end-1));

% No change of sign at the beginning
sigdif2 = [0; sigdif2];
sigdifFramed2 = buffer(sigdif2, winLen, winOverlap,'nodelay');
sigdifWindowed2 = diag(sparse(wHamm)) * abs(sigdifFramed2);
zcr2 = sum(sigdifWindowed2, 1)/(2*winLen);

t = (0:length(speechSignal2)-1)/Fs;
delay = (winLen - 1)/2;
plot(t(delay+1:end-delay), zcr2/max(zcr2), 'g');
legend({'Music cut','Short-Time Energy','Zero Crossing Rate'});

hold off;