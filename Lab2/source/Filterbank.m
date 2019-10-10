function [h,g] = Filterbank(M)

L=2*M;

h = zeros(L,M);
for n=0:L-1
    for k=0:M-1
        h(n+1,k+1) = sin((n+0.5)*pi/(2*M))*sqrt(2/M)*cos((2*n+M+1)*(2*k+1)*pi/(4*M));
    end
end

g = zeros(L,M);
for n=1:L
    for k=1:M
        g(n,k) = h(2*M+1-n,k);
    end
end


end