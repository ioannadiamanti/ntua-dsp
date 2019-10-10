function [Ttm,Tnm] = IndividualMaskingThresholds(Ptm,Pnm,b)
N = length(Ptm);    
Ttm = zeros(N,N);
Tnm = zeros(N,N);
for i=1:N
    for j=1:N
        Db = b(i)-b(j);
        if Db>=-3 && Db<-1
            if Ptm(j)>0
                SFtm = 17*Db-0.4*Ptm(j)+11;
                Ttm(i,j)= Ptm(j)-0.275*b(j) + SFtm - 6.025;
            end
            if Pnm(j)>0
                SFnm = 17*Db-0.4*Pnm(j)+11;
                Tnm(i,j)=  Pnm(j) - 0.175*b(j)+SFnm-2.025;
            end
        elseif Db>=-1 && Db<0
            if Ptm(j)>0
                SFtm = (0.4*Ptm(j)+6)*Db;
                Ttm(i,j)= Ptm(j)-0.275*b(j) + SFtm - 6.025;
            end
            if Pnm(j)>0
                SFnm = (0.4*Ptm(j)+6)*Db;
                Tnm(i,j)=  Pnm(j) - 0.175*b(j)+SFnm-2.025;
            end
        elseif Db>=0 && Db<1
            if Ptm(j)>0
                SFtm = -17*Db;
                Ttm(i,j)= Ptm(j)-0.275*b(j) + SFtm - 6.025;
            end
            if Pnm(j)>0
                SFnm = -17*Db;
                Tnm(i,j)=  Pnm(j) - 0.175*b(j)+SFnm-2.025;
            end
        elseif Db>=1 && Db<8
            if Ptm(j)>0
                SFtm = (0.15*Ptm(j)-17)*Db-0.15*Ptm(j);
                Ttm(i,j)= Ptm(j)-0.275*b(j) + SFtm - 6.025;
            end
            if Pnm(j)>0
                SFnm = (0.15*Pnm(j)-17)*Db-0.15*Pnm(j);
                Tnm(i,j)=  Pnm(j) - 0.175*b(j)+SFnm-2.025;
            end
        end 
    end
end
end