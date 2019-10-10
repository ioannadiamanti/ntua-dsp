function Tg = GlobalMaskingThreshold(Ttm,Tnm,Tq)
    Tg = 10*log10(10.^(0.1*Tq)+sum(10.^(0.1*Ttm),2) + sum(10.^(0.1*Tnm),2));
end