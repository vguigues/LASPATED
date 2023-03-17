



%Districts
lambda103D=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_District.txt');
LHD=lambda103D;
LDP=zeros(48*7,1);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:160
            for p=1:3
                thisest=thisest+LHD((g-1)*48+t,(p-1)*160+r);
            end
        end
        LDP((g-1)*48+t)=thisest;
    end
end


plot(estTarget,'-k');
%Districts
plot(LDP,'--r');