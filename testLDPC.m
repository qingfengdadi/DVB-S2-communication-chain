clear; close all;
H = [0 1 0 1 1 0 0 1;
     1 1 1 0 0 1 0 0;
     0 0 1 0 0 1 1 1;
     1 0 0 1 1 0 1 0];

tic;
graph = buildTannerGraph(H);
msg_tx = [1 0 0 1 0 1 0 1];
msg_rx = [1 1 0 1 0 1 0 1];

g2 = LdpcSoftDecoder(msg_rx, H, graph, 1, 5);
toc


%%
N0 = 1;
symb_rx = msg_tx;
Lq = -2*symb_rx'/N0;

bit_out = ones(1,length(symb_rx));
sizeH = size(H);

C = struct('f', 0 ,'bit', symb_rx(1)');
it = 1;

for i = 2:sizeH(2)
    C(end+1) = struct('f', 0, 'bit', symb_rx(i)');
end
C_save = C;
%init
F = struct('c', find(H(1,:)),'bit', Lq(find(H(1,:))));
for i = 2:sizeH(1)
    temp = find(H(i,:));
    F(end+1) = struct('c', temp, 'bit', Lq(temp));
end

while(sum(mod(bit_out*H',2)) ~= 0 && it < 10)


    if (it >1)
        for i = 1:sizeH(1)
            F(i).bit = [];
        end
        
        for i = 1:sizeH(2)
            f_ = C(i).f;
            bit_ = C(i).bit;
            for k = 2:length(f_)
                hihi = sum([bit_(1:k-1) bit_(k+1:end)]);
                F(f_(k)).bit = [F(f_(k)).bit hihi];
            end
        end
        C = C_save;
    end
    
    for j = 1:sizeH(1)
        c_ = F(j).c;
        bit_ = F(j).bit;
        for k = 1:length(c_);
            C(c_(k)).f = [C(c_(k)).f j];
            product = prod(sign([bit_(1:k-1) bit_(k+1:end)]));
            alpha = abs([bit_(1:k-1) bit_(k+1:end)]);
            C(c_(k)).bit = [C(c_(k)).bit product*min(alpha)];
        end
    end

    for i = 1:sizeH(2)
        c_ = C(i).bit;
        if (sum(c_))<0
            bit_out(i) = 1;
        else
            bit_out(i) = 0;
        end        
    end
    it = it + 1;
end

bit_hihi = bit_out(length(bit_out)/2+1:end)';
it
bit_tx = symb_rx;
err = sum(bit_tx-bit_hihi)