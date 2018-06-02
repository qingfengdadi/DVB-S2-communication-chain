clear; close all
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));

H = [1 1 0 1 1 0 0 1 0 0; 
     0 1 1 0 1 1 1 0 0 0;
     0 0 0 1 0 0 0 1 1 1;
     1 1 0 0 0 1 1 0 1 0;
     0 0 1 0 0 1 0 1 0 1];
 
% Build H such that H = [I P']
H(2,:) = xor(H(2,:),H(5,:));
H(1,:) = xor(H(1,:),H(3,:));
H(1,:) = xor(H(1,:),H(2,:));
temp = H(4,:);
H(4,:) = H(3,:);
H(3,:) = H(5,:);
H(5,:) = temp;
H(5,:) = xor(H(5,:), H(1,:));
H(2,:) = xor(H(2,:), H(5,:));
temp = H(2,:);
H(2,:) = H(5,:);
H(5,:) = temp;

% Build G from H
P = H(:,6:end)';
G = [P eye(5)];

message_vector = [0 1 1 0 1];
u = mod(message_vector*G,2) % code vector = [0 1 0 0 0 0 1 1 0 1]
r = [0 1 0 0 0 0 1 1 1 1]; % received vector

graph = buildTannerGraph(H);
u_ = LdpcHardDecoder(r, H, graph, 15) % corrected vector
