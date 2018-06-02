function [n, df] = cfoEstimate(rx, p, T, K)
%{
    input
    rx: received symbols
    p:  pilot symbols
    T:  symbol duration
    K:  K-factor

    output
    n:  arrival time estimate
    df: cfo estimate
%}

N = length(p);
L = length(rx);
D = zeros(K, L-N+1);

for k = 1 : K
    for n = 1 : L - N + 1
        tmp = 0;
        for l = k : N-1
            tmp1 = conj(rx(n+l)) * p(l+1);
            
            tmp2 = conj(rx(n+l-k)) * p(l-k+1);
            
            tmp = tmp + tmp1 * conj(tmp2);
        end
        D(k,n) = tmp;
        
    end 
    D(k,:) = D(k,:)/(N-k);
end

% Time of arrival estimate
tmp = sum(abs(D),1);
[~, n] = max(tmp);

% CFO estimate
% df = 0;
% for k = 1 : K
%     df = df + angle(D(k,n))/(2*pi*k*T);
% end
% df = (-1/K) * df

k = (1:K).' ;
df = sum(angle(D(k,n))./(2*pi*k*T), 1);
df = - df/K;

end