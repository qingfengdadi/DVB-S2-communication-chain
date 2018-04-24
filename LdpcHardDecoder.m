function [graph] = LdpcHardDecoder(bits_rx, graph)
N = length(bits_rx);
% first step: initialize the v_nodes
for i=1:N
    graph.v_nodes(i).msg = bits_rx(i);
end
% then communicate the msg to c_nodes
for i=1:N
    c_nodes = [graph.v_nodes(i).c_nodes];
    for j=1:length(c_nodes)
        graph.c_nodes(c_nodes(j)).msg = graph.v_nodes(i).msg;
    end
end
end

