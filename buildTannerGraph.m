function [graph] = buildTannerGraph(H)
% example of use for a c_node:
% to access the n_th c_node:    c = graph.c_nodes(n)
% each c_node has two attributes:
% 1- the msg sent by that c_node : c.msg
% 2- the v_nodes connected to that c_node: v = c.v_nodes
% to access a specific v_node: choose the m_th relative connection with the c_node c:
% vm = c.v_nodes(m), which is a structure containing two information:
% 1- the general index number of the m_th connected v_node:    vm.num
% 2- the msg received by that connected v_node: vm.msg

[m,n] = size(H);
for i=1:m
    c_to_v_nodes = find(H(i,:));
    for j=1:length(c_to_v_nodes)
        c_nodes(i).msg = 0; % msg sent by c_node
        c_nodes(i).v_nodes(j).num = c_to_v_nodes(j);
        c_nodes(i).v_nodes(j).msg = 0; % msg received by c_node
    end
end
    
for i=1:n
    v_to_c_nodes = find(H(:,i));
    for j=1:length(v_to_c_nodes)
        v_nodes(i).msg = 0; % msg sent by v_node
        v_nodes(i).c_nodes(j).num = v_to_c_nodes(j);
        v_nodes(i).c_nodes(j).msg = 0; % msg received by v_node
    end
end

graph = struct();
graph.c_nodes = c_nodes;
graph.v_nodes = v_nodes;
end


