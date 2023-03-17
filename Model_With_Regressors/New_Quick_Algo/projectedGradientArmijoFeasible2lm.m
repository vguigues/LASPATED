
function [beta,fVal]=projectedGradientArmijoFeasible2lm(nbObservations,nbCalls,beta,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBetaA,Groups,durations,sizex,whichgroup,lbounds,gnb,cnb,nbC)

k=1;
fVal=[];
bparam=2;
betak=bparam;

while (k<=iterMax)
      k
      [fold]=oracleObjectiveModel2lm(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBetaA,gnb,cnb);
      [gradient]=oracleGradientModel2lm(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,sizex,Groups,whichgroup,weight,durations,indexBetaA,gnb,cnb);                                
      [z]=projectionRegressorlm(regressor,indexBetaA,nbLandTypes,C,D,T,R,epsilon,sizex,x-betak*gradient,lbounds,gnb,cnb,nbC,Groups);
      bool=1;
      j=0;
      rhs=gradient'(x-z); 
      while (bool==1)
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel2lm(nbObservations,nbCalls,zAux,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBetaA,gnb,cnb);
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      k=k+1;
      betak=bparam/(2^j);
end

