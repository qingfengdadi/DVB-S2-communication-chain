addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));

clear all; close all;
% Generate a stream tx of random bits
N = 1000; % length of the stream
fs = 10^6;
Nbps = 10;
modulation = 'qam'; % or pam
bit_tx = randi(2,N,1)-1;
symb_tx = mapping(bit_tx,Nbps,modulation);

Fsampling = 8e5;
beta = 0.5;
Tsymb = 1e-6;

f = -Fsampling/2:1000:Fsampling/2;
HRRC = HRRC(f,Tsymb,beta);
plot(f,HRRC)