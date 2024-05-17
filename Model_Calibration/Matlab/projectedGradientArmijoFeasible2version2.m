

function [x,fVal]=projectedGradientArmijoFeasible2version2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,iterMax,sigma,epsilon,indexBeta,durations,sizex)

k=1;
fVal=[];
bparam=2;
betak=bparam;
lbounds=zeros(sizex,1);
%upper=10^(50);
%lower=-10^(50)
[x]=projectionRegressor(regressor,indexBeta,nbRegressors,C,D,T,R,epsilon,sizex,x,lbounds);


while (k<=iterMax)
      [fold]=oracleObjectiveModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,durations,indexBeta);
      fold
      [gradient]=oracleGradientModel2version2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,sizex,durations,indexBeta);
      betak
      sum(gradient)
      disp('sum x-z')
      sum(x-z)
      x-betak*gradient
      [z]=projectionRegressor(regressor,indexBeta,nbRegressors,C,D,T,R,epsilon,sizex,x-betak*gradient,lbounds);
      
      bool=1;
      j=0;
      if (k>=2)
          gradient

      end
      aux1=[];
      aux2=[];
      for i=1:sizex
          x(i)-z(i)
          if (gradient(i)*(x(i)-z(i))>=0)
              aux1=[aux1;gradient(i)*(x(i)-z(i))];
          else
              aux2=[aux2;-gradient(i)*(x(i)-z(i))];
          end
      end
      [sortedaux1]=sort(aux1);
      [sortedaux2]=sort(aux2);
      rhs=0;
      for i=1:length(aux1)
          rhs=rhs+sortedaux1(i);
      end
      rhs2=0;
      for i=1:length(aux2)
          rhs2=rhs2+sortedaux2(i);
      end
      rhs=rhs-rhs2
      while (bool==1)
          %bool
          %j
          %rhs
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel2(nbObservations,nbArrivals,zAux,C,D,T,R,regressor,nbRegressors,durations,indexBeta);        
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      k=k+1;
      %upper=min(upper,f);
      betak=bparam/(2^j);
end
