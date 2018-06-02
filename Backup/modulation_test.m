%% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));
clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RANDOM BITS
N = 1000; % bit stream length
bits_tx = randi(2,N,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 8; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols
Nbps = 4; % number of bits per symbol
modulation = 'qam'; % type of modulation 
beta = 0.3; % roll-off factor

%% mapping
signal = mapping(bits_tx,Nbps,modulation);
% figure;
% plot(signal,'*');
%% upsampling
signal_tx = upsample(signal,M);

%% implementation of transfer function (freq)
RRCtaps = 365;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrc = HRC(f, Tsymb, beta);
Hrrc = sqrt(Hrc);
% figure;
% semilogy(f,Hrrc);

%% implementation of transfer function (time)
deltat = 1/fsampling;
t = (-(RRCtaps-1)/2:(RRCtaps-1)/2)*deltat;
g = ifftshift(ifft(fftshift(Hrrc)));
h = conv(g, g);
h_norm = h/max(h);
h_time = h_norm(RRCtaps:end-RRCtaps+1);
figure;
semilogy(t, h_time);

%% Symbol interference illustration
figure;
plot(t,h(RRCtaps:end-RRCtaps+1));
hold on
plot(t+Tsymb,h(RRCtaps:end-RRCtaps+1));
hold on
plot(t-Tsymb,h(RRCtaps:end-RRCtaps+1));
hold on
plot(t+2*Tsymb,h(RRCtaps:end-RRCtaps+1));
hold on
plot(t-2*Tsymb,h(RRCtaps:end-RRCtaps+1));
grid on;

%%
