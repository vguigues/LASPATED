

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the gradient of the objective for calibrating the model for 
%the process of calls. Model 1 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient: gradient of the objective function at lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gradient]=oracleGradientModel2(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,sizex,Groups,whichgroup,weight,durations,indexBeta,gnb,cnb)

gradient=zeros(sizex,1);
for m=1:length(Groups{1,gnb})
    thisday=Groups{1,gnb}(m).d;
    thistime=Groups{1,gnb}(m).t;
    for i=1:R
        rates=0;
        for j=1:(1+nbLandTypes)
            rates=rates+x(indexBeta(cnb,thisday,thistime,j))*regressor(j,i);
        end
        for j=1:(1+nbLandTypes)
            gradient(indexBeta(c,d,t,j))=gradient(indexBeta(c,d,t,j))+nbObservations(cnb,thisday,thistime,i)*regressor(j,i)-(nbCalls(cnb,thisday,thistime,i)*regressor(j,i)/rates);
        end
    end
end

gnb,cnb

if (length(Groups)>0)
    for index1=1:length(Groups{1,gnb})
        d1=Groups{1,gnb}(index1).d;
        t1=Groups{1,gnb}(index1).t;
        for index2=1:length(Groups{1,gnb})
            d2=Groups{1,gnb}(index2).d;
            t2=Groups{1,gnb}(index2).t;    
            if (index2~=index1)
                for j=1:(1+nbLandTypes)
                    gradient(indexBeta(cnb,d1,t1,j))=gradient(indexBeta(cnb,d2,t2,j))+(2*weight/durations(t1))*((gradient(indexBeta(cnb,d1,t1,j))/durations(t1))-(gradient(indexBeta(cnb,d2,t2,j))/durations(t2)));
                end
            end
        end
    end
end



