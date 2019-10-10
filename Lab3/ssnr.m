function result = ssnr(sigWindowed,noiseWindowed)
temp = log10(sum(sigWindowed.^2,1)./sum(noiseWindowed.^2));
index35 = find(temp>35);
index0 = find(temp<-10);
temp(index35) = 35;
temp(index0) = 0;
result = 10/(size(temp,2)-length(index0)) * sum(temp);
end