function Ptm = findToneMaskers(p)
St = zeros(length(p),1);
for k=1:length(p)
    if k < 3 || k>250
        St(k) = 0;
    elseif k<63
        St(k) = (p(k)>p(k-1) & p(k)>p(k+1) & p(k)>p(k+2)+7 & p(k)>p(k-2)+7);
    elseif k<127
        St(k) = (p(k)>p(k-1) & p(k)>p(k+1) & p(k)>p(k+2)+7 & p(k)>p(k-2)+7 & p(k)>p(k+3)+7 & p(k)>p(k-3)+7);
    else
        St(k) = (p(k)>p(k-1) & p(k)>p(k+1) & p(k)>p(k+2)+7 & p(k)>p(k-2)+7 & p(k)>p(k+3)+7 & p(k)>p(k-3)+7 & p(k)>p(k+4)+7 & p(k)>p(k-4)+7 & p(k)>p(k+5)+7 & p(k)>p(k-5)+7 & p(k)>p(k+6)+7 & p(k)>p(k-6)+7);
    end
end 

Ptm = zeros(length(p),1);
for k = 1:length(p)
    if St(k)==1
        Ptm(k) = 10*log10(10^(0.1*p(k-1))+10^(0.1*p(k))+10^(0.1*p(k+1)));
    end
end
end