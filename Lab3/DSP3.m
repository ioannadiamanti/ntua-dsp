clc;
clear all;
close all;

p = which('DSP3');
idcs = strfind(p,'\');
p = p(1:idcs(end)-1);
cd(p);
addpath(genpath('Material\LogarithmicPolarPlot'));

%%
%Μέρος 1ο
%1.4
theta_s = pi/2;
f = 2000;
w = 2*pi*f;
theta = linspace(0,pi,f);
moires=(180/pi)*theta;
c = 340;

%d = 8cm, N = 4,8,12,16
d=0.08;
N=4;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
plot(moires,B);
hold on;

N=8;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);

N=12;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);


N=16;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);

ylabel('log(|b=B(ω,θ)|)');
xlabel('theta (degrees)');
title('d = 8cm, N = 4,8,12,16, θs = 90');
legend('N=4','N=8','N=12','N=16');
set(gca, 'YScale', 'log');
hold off;

%N = 8, d = 8,12,16,20 cm
N=8;
d=0.08;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
plot(moires,B);
hold on;

d=0.12;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);

d=0.16;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);

d=0.2;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(moires,B);

ylabel('log(|b=B(ω,θ)|)');
xlabel('theta (degrees)');
title('N = 8, d = 8,12,16,20 cm, θs = 90');
legend('d=8cm','d=12cm','d=16cm','d=20cm');
set(gca, 'YScale', 'log');
hold off;

%N=8, d = 8cm, theta_s = 90 , theta = 0,45,60 degrees
N=8;
d=0.08;
f = linspace(0,8000,1000);
w = 2*pi*f;

theta = 0;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
plot(f,B);
hold on;

theta = pi/4;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(f,B);

theta = pi/3;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
plot(f,B);

ylabel('log(|b=B(ω,θ)|)');
xlabel('Frequency (Hz)');
title('N = 8, d = 8cm, θs = 90, θ = 0, 45, 60');
legend('θ= 0 ','θ=45','θ=60');
set(gca, 'YScale', 'log');
hold off;

%N=8, d=0.08, f = 2kHz, θs = 0,45,90
N=8;
d=0.08;
f=2000;
w=2*pi*f;
theta = linspace(-pi,pi,f);

theta_s = 0;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
semilogr_polar(theta,B);
title('Polar Diagram of |log(B(ω,θ))| for θs = 0');

theta_s = pi/4;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
semilogr_polar(theta,B);
title('Polar Diagram of |log(B(ω,θ))| for θs = 45');

theta_s = pi/2;
B = abs((1/N)*sin((N/2)*(w/c)*d*(cos(theta)-cos(theta_s)))./sin(0.5*(w/c)*d*(cos(theta)-cos(theta_s))));
figure;
semilogr_polar(theta,B);
title('Polar Diagram of |log(B(ω,θ))| for θs = 90');

%%
%Μέρος 2ο
%2.1 Beamforming σε προσομοιωμένα σήματα

addpath(genpath('Material\MicArrayRealSignals'));
addpath(genpath('Material\MicArraySimulatedSignals'));


N=7;
d=0.08;
theta_voice = pi/4;
theta_noise = 3*pi/4;
f_noise = linspace(500,2500,1000);

[source,Fs] = audioread('Material\MicArraySimulatedSignals\source.wav');
sensor0 = audioread('Material\MicArraySimulatedSignals\sensor_0.wav');
sensor1 = audioread('Material\MicArraySimulatedSignals\sensor_1.wav');
sensor2 = audioread('Material\MicArraySimulatedSignals\sensor_2.wav');
sensor3 = audioread('Material\MicArraySimulatedSignals\sensor_3.wav');
sensor4 = audioread('Material\MicArraySimulatedSignals\sensor_4.wav');
sensor5 = audioread('Material\MicArraySimulatedSignals\sensor_5.wav');
sensor6 = audioread('Material\MicArraySimulatedSignals\sensor_6.wav');


%A) Delay-and-sum beamforming
n= (0:6);
pn = (n - (N-1)/2)*d;
w = linspace(0, 2*pi, length(source))*Fs;
as = cos(theta_voice);
ks = w/c*as;
dks = zeros(N,length(source));
for i=1:N
    dks(i,:) = exp(-1j*ks*pn(i));
end
H = 1/N * dks.';
f = [sensor0 sensor1 sensor2 sensor3 sensor4 sensor5 sensor6];
F = fft(f);
Y = sum(H.*F,2);
y = ifft(Y);
y = real(y);

%Plots

%Plot source signal
figure;
t = (0:length(source)-1)/Fs;
plot(t,source);
title('Source Signal');
ylabel('Amplitude');
xlabel('Time(sec)');

%Plot central microphone signal
figure;
t = (0:length(sensor3)-1)/Fs;
plot(t,sensor3);
title('Central Signal');
ylabel('Amplitude');
xlabel('Time(sec)');

%Plot Beamformer Output
figure;
t = (0:length(y)-1)/Fs;
plot(t,y);
title('Beamformer Output Signal');
ylabel('Amplitude');
xlabel('Time(sec)');

%Spectograms
L = Fs*0.03;
P = 2/3 * L;


figure;
spectrogram(source,L,P,L,Fs,'yaxis');
title('Spectrogram of Source Signal');

figure;
spectrogram(sensor3,L,P,L,Fs,'yaxis');
title('Spectrogram of Central Microphone Signal');

figure;
spectrogram(y,L,P,L,Fs,'yaxis');
title('Spectrogram of Beamformer Output Signal');

noise3 = source - sensor3;
noisey = source - y;

SNRs3 = snr(source,noise3);
SNRy = snr(source,noisey);

audiowrite('sim_ds.wav',y,Fs);

%B) Single-Channel Wiener Filtering
from = 0.47;
to = 0.5;
x = sensor3(from*Fs:to*Fs);
s = source(from*Fs:to*Fs);
v = x-s;

%10ms window
win = Fs*0.01;
winoverlap = win/2;
NFFT = Fs*0.03+1;

[px,f] = pwelch(x,win,winoverlap,NFFT,Fs,'twosided');
[pv,~] = pwelch(v,win,winoverlap,NFFT,Fs,'twosided');

Hw = 1 - pv./px;

figure;
semilogy(f,abs(Hw));
xlim([0 8000]);
ylabel('Amplitude (dB)');
xlabel('Frequency (Hz)');
title('Wiener Filter Response');

nsd = abs(1-Hw).^2;
figure;
semilogy(f,nsd);
xlim([0 8000]);
ylabel('Amplitude (dB)');
xlabel('Frequency (Hz)');
title('Speech Distortion Index');

Yw = Hw.*fft(x);
yw = ifft(Yw);

ps = pwelch(s,win,winoverlap,NFFT,Fs,'twosided');
pyw = pwelch(yw,win,winoverlap,NFFT,Fs,'twosided');
figure;
semilogy(f,ps);
hold on;
semilogy(f,px);
semilogy(f,pyw);
semilogy(f,pv);
xlim([0 8000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Power Spectrum');
legend('Power Spectrum of s(t)','Power Spectrum of x(t)','Power Spectrum of yw(t)','Power Spectrum of v(t)');
hold off;

SNRx = snr(x,v);
SNRyw = snr(s,real(s-yw));
SNRbeamformer = snr(s,s-y(from*Fs:to*Fs));

py = pwelch(y(from*Fs:to*Fs),win,winoverlap,NFFT,Fs,'twosided');

figure;
semilogy(f,ps);
hold on;
semilogy(f,px);
semilogy(f,pyw);
semilogy(f,py);
xlim([0 8000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Power Spectrum');
legend('Power Spectrum of s(t)','Power Spectrum of x(t)','Power Spectrum of yw(t) (wiener)','Power Spectrum of y(t) (beamformer)');
hold off;

%2.2 Beamforming σε πραγματικά σήματα
N = 7;
d = 0.04;
theta = pi/4;

[source,Fs] = audioread('Material\MicArrayRealSignals\source.wav');
sensor0 = audioread('Material\MicArrayRealSignals\sensor_0.wav');
sensor1 = audioread('Material\MicArrayRealSignals\sensor_1.wav');
sensor2 = audioread('Material\MicArrayRealSignals\sensor_2.wav');
sensor3 = audioread('Material\MicArrayRealSignals\sensor_3.wav');
sensor4 = audioread('Material\MicArrayRealSignals\sensor_4.wav');
sensor5 = audioread('Material\MicArrayRealSignals\sensor_5.wav');
sensor6 = audioread('Material\MicArrayRealSignals\sensor_6.wav');

%A) Delay-and-sum beamforming
n= (0:6);
pn = (n - (N-1)/2)*d;
w = linspace(0, 2*pi, length(source))*Fs;
as = cos(theta);
ks = w/c*as;
dks = zeros(N,length(source));
for i=1:N
    dks(i,:) = exp(-1j*ks*pn(i));
end
H = 1/N * dks.';
f = [sensor0 sensor1 sensor2 sensor3 sensor4 sensor5 sensor6];
F = fft(f);
Y2 = sum(H.*F,2);
y2 = ifft(Y2);
y2 = real(y2);

audiowrite('real_ds.wav',y2,Fs);

%Plot source signal
L = Fs*0.03;
P = 2/3*L;

figure;
t = (0:length(source)-1)/Fs;
plot(t,source);
xlabel('Time(sec)');
ylabel('Amplitude');
title('Source Signal');


figure;
spectrogram(source,L,P,L,Fs,'yaxis');
title('Spectrogram of Source Signal');

figure;
t = (0:length(sensor3)-1)/Fs;
plot(t,sensor3);
xlabel('Time(sec)');
ylabel('Amplitude');
title('Central Microphone Signal');

figure;
spectrogram(sensor3,L,P,L,Fs,'yaxis');
title('Spectrogram of Central Microphone Signal');

figure;
t = (0:length(y2)-1)/Fs;
plot(t,y2);
xlabel('Time(sec)');
ylabel('Amplitude');
title('Beamformer Output Signal');


figure;
spectrogram(y2,L,P,L,Fs,'yaxis');
title('Spectrogram of Beamformer Output Signal');

noise = source-sensor3;

% Windowing the signals
sigWindowed = buffer(sensor3,100);

%WSS noise - taking noise only from one frame
noiseWindowed = noise(1:100);

ssnr_sensor3 = ssnr(sigWindowed,noiseWindowed);

noise = source-y2;

figure;
t = (0:length(noise)-1)/Fs;
plot(t,noise);
title('Diffuse Noise Field');
xlabel('Time(Sec)');
ylabel('Amplitude');

% Windowing the signals
sigWindowed = buffer(y2,100);

%WSS noise - taking noise only from one frame
noiseWindowed = noise(1:100);

ssnr_beamformer = ssnr(sigWindowed,noiseWindowed);

%B) Post - filtering using Wiener Filter
winlen = Fs*0.03;
win = hamming(winlen);
winOverlap = 2*winlen/3;
noise = sensor3(1:winlen);

pnoise = pwelch(noise,winlen,winOverlap,winlen,Fs,'twosided');

yw2 = post_filter(y2,pnoise,winlen,winlen/3,winlen/6,winlen,Fs);

audiowrite('real_mmse.wav',yw2,Fs);



%plots

figure;
t = (0:length(source)-1)/Fs;
plot(t,source);
title('Source signal');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
t = (0:length(sensor3)-1)/Fs;
plot(t,sensor3);
title('Central Microphone Signal');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
t = (0:length(y2)-1)/Fs;
plot(t,y2);
title('Wiener Input signal');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
t = (0:length(yw2)-1)/Fs;
plot(t,yw2);
title('Wiener Output signal');
xlabel('Time(sec)');
ylabel('Amplitude');

%spectrograms

L = Fs*0.03;
P = 2/3 * L;


figure;
spectrogram(source,L,P,L,Fs,'yaxis');
title('Spectrogram of Source Signal');

figure;
spectrogram(sensor3,L,P,L,Fs,'yaxis');
title('Spectrogram of Central Microphone Signal');

figure;
spectrogram(y2,L,P,L,Fs,'yaxis');
title('Spectrogram of Wiener Input Signal');


figure;
spectrogram(yw2,L,P,L,Fs,'yaxis');
title('Spectrogram of Wiener Output Signal');

noise = source - y2;
noiseWindowed = noise(1:100);
sigWindowed = buffer(y2,100);

ssnr_input_wiener = ssnr(sigWindowed,noiseWindowed);

noise = source - yw2.';
noiseWindowed = noise(1:100);
sigWindowed = buffer(yw2,100);
ssnr_output_wiener = ssnr(sigWindowed,noiseWindowed);

sum_ssnr = 0;

for i=1:size(f,2)
    noise = source - f(:,i);
    noiseWindowed = noise(1:100);
    sigWindowed = buffer(f(:,i),100);
    sum_ssnr = sum_ssnr + ssnr(sigWindowed,noiseWindowed);
end
meanssnr = sum_ssnr/size(f,2);
