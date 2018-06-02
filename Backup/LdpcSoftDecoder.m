function [decoded_bits] = LdpcSoftDecoder(bits_rx, H, graph, noiseVar, numIt)
[M,N] = size(H);
it = 0;
decoded_bits = bits_rx;
Lci = (-2/noiseVar)*bits_rx;
% first step: 
% initialize the v_nodes
% for i=1:N
%     graph.v_nodes(i).msg = Lci;
% end
% then communicate the msg to c_nodes
for i=1:N
    c_nodes = graph.v_nodes(i).c_nodes; % c_nodes connected to a given v_node
    for j=1:length(c_nodes)
        idx = c_nodes(j).num; % get the global c_node index
        nums = [graph.c_nodes(idx).v_nodes(:).num];  % v_nodes connected the c_node
        for k=1:length(nums)
            % to get relative index of v_node with respect to c_node
            if nums(k) == i
                graph.c_nodes(idx).v_nodes(k).msg = Lci(i);
            end
        end
    end
end

while ((it < numIt) && (length(find(mod(decoded_bits*H',2))) ~= 0))
    % second step:
    % calculate the response of the c_nodes using Gallager formula 
    % (in the log domain)
    for i=1:M
        v_nodes = graph.c_nodes(i).v_nodes; % v_nodes connected to a given c_node
        for j=1:length(v_nodes)
            v_nodes_no_j = [v_nodes(1:j-1) v_nodes(j+1:end)];
            alpha_no_j = min(abs([v_nodes_no_j(:).msg]));
            prodXi = prod(sign([v_nodes_no_j(:).msg]));
            idx = v_nodes(j).num;
            nums = [graph.v_nodes(idx).c_nodes(:).num]; % c_nodes connected the v_node
            Lrji = prodXi * alpha_no_j;
            for k=1:length(nums)
                % to get relative index of v_node with respect to c_node
                if nums(k) == i
                    graph.v_nodes(idx).c_nodes(k).msg = Lrji;
                end
            end
        end
    end

    % third step:
    % v_nodes take a decision based on the response from c_node and original
    % msg using a majoroty vote
    for i=1:N
        c_nodes = graph.v_nodes(i).c_nodes;
        softDecision = Lci(i) + sum([c_nodes(:).msg]);
        for j=1:length(c_nodes)
            idx = c_nodes(j).num; % get the global c_node index
            nums = [graph.c_nodes(idx).v_nodes(:).num];  % v_nodes connected the c_node
            Lqij = softDecision - c_nodes(j).msg;
            for k=1:length(nums)
                % to get relative index of v_node with respect to c_node
                if nums(k) == i
                    graph.c_nodes(idx).v_nodes(k).msg = Lqij;
                end
            end
        end
        if softDecision < 0
            decoded_bits(i) = 1;
        else
            decoded_bits(i) = 0;
        end
    end
    it = it + 1;
end

end