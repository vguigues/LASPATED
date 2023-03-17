
function [fold]=oracleObjectiveModel2lm(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBeta,gnb,cnb)

gnb,cnb
f=0;
for m=1:length(Groups{1,gnb})
    thisday=Groups{1,gnb}(m).d;
    thistime=Groups{1,gnb}(m).t;
    for i=1:R
        rates=0;
        for j=1:(1+nbLandTypes)
            rates=rates+x(indexBeta(cnb,thiday,thistime,j))*regressor(j,i);
        end
        f=f+nbObservations(cnb,thiday,thistime,i)*rates-nbCalls(cnb,thiday,thistime,i)*log(rates);
    end
end

if (length(Groups)>0)
    for index1=1:length(Groups{1,gnb})
        d1=Groups{1,gnb}(index1).d;
        t1=Groups{1,gnb}(index1).t;
        for index2=1:length(Groups{1,gnb})
            if (index2~=index1)
                d2=Groups{1,gnb}(index2).d;
                t2=Groups{1,gnb}(index2).t;
                for j=1:(1+nbLandTypes)
                    f=f+(weight/2)*((x(indexBeta(cnb,d1,t1,j))/durations(t1))-(x(indexBeta(cnb,d2,t2,j))/durations(t2)))^2;
                end
            end
        end
    end
end
