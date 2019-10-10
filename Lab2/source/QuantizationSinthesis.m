function [output,bits] = QuantizationSinthesis(x,h,g,Tg,M,adaptive)

R = 2^16;
bits = 0;
for k=1:M
    v = conv(h(:,k),x);
    y = downsample(v,M);
    if adaptive == 1
        Bk = ceil(log2(R/min(Tg(8*(k-1)+1:8*k)))-1);
        ymax = max(y);
        ymin = min(y);
    else
        Bk = 8;
        ymax=1;
        ymin=-1;
    end
    bits = bits + Bk*length(y);
    Delta = abs(ymax-ymin)/(2^Bk);
    y_hat = zeros(length(y),1);
    for i=1:length(y)
        temp = floor((y(i)-ymin)/Delta);
        if temp == 2^Bk
            temp = temp-1;
        end
        y_hat(i) = temp*Delta+Delta/2+ymin;
    end
    
    w(:,k) = upsample(y_hat,M);
    
end
    %bits = bits/M;
    output = conv(g(:,1),w(:,1));
for k=2:M
    output = output + conv(g(:,k),w(:,k));
end
end