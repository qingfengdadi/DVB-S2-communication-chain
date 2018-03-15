addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));

% Generate a stream tx of random bits
N = 1000; % length of the stream
Nbps = 10;
modulation = 'qam'; % or pam
bit_tx = randi(2,N,1)-1;
symb_tx = mapping(bit_tx,Nbps,modulation);

