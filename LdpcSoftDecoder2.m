function [decoded_bits, it] = LdpcSoftDecoder2(symb_rx, H, N0, numIt)
Lq = -2*symb_rx'/N0;

decoded_bits = ones(1,length(symb_rx));
[m,n] = size(H);

it = 1;
for i = 1:n
    C(i) = struct('f', 0, 'bit', symb_rx(i)');
end
C_save = C;

for i = 1:m
    indexes = find(H(i,:));
    F(i) = struct('c', indexes,'bit', Lq(indexes));
end     

while(sum(mod(decoded_bits*H',2)) ~= 0 && it < numIt)
    if (it >1)
        for i = 1:m
            F(i).bit = [];
        end
        
        for i = 1:n
            f_ = C(i).f;
            bit_ = C(i).bit;
            for k = 2:length(f_)
                res = sum([bit_(1:k-1) bit_(k+1:end)]);
                F(f_(k)).bit = [F(f_(k)).bit res];
            end
        end
        C = C_save;
    end
    
    for j = 1:m
        c_ = F(j).c;
        bit_ = F(j).bit;
        for k = 1:length(c_)
            C(c_(k)).f = [C(c_(k)).f j];
            product = prod(sign([bit_(1:k-1) bit_(k+1:end)]));
            alpha = abs([bit_(1:k-1) bit_(k+1:end)]);
            C(c_(k)).bit = [C(c_(k)).bit product*min(alpha)];
        end
    end

    for i = 1:n
        c_ = C(i).bit;
        if (sum(c_))<0
            decoded_bits(i) = 1;
        else
            decoded_bits(i) = 0;
        end        
    end
    it = it + 1;
end

decoded_bits = decoded_bits(length(decoded_bits)/2+1:end)';