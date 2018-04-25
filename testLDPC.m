
H = [0 1 0 1 1 0 0 1;
     1 1 1 0 0 1 0 0;
     0 0 1 0 0 1 1 1;
     1 0 0 1 1 0 1 0];

tic;
graph = buildTannerGraph(H);
msg_tx = [1 0 0 1 0 1 0 1];
msg_rx = [1 1 0 1 0 1 0 1];

decoded_msg = LdpcHardDecoder(msg_rx, H, graph, 2)
tocms