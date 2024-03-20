

function [x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,iterMax,sigma,epsilon,indexBeta,durations,sizex)

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
      [gradient]=oracleGradientModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,sizex,durations,indexBeta);
      [z]=projectionRegressor(regressor,indexBeta,nbRegressors,C,D,T,R,epsilon,sizex,x-betak*gradient,lbounds);
      bool=1;
      j=0;
      % diff_aux = x-z;
      % for c=1:C
      %     for d=1:D
      %         for t=1:T
      %             for j=1:nbRegressors
      %                 fprintf("c%i d%i t%i j%i | x = %.8f grad = %.8f x-z = %.8f\n", c-1,d-1,t-1,j-1, x(indexBeta(c,d,t,j)), gradient(indexBeta(c,d,t,j)), diff_aux(indexBeta(c,d,t,j)));
      %             end
      %         end
      %     end
      % end
      % fprintf("rhs = %.15f\n",gradient'*(x-z));
      % input("");
      rhs=gradient'*(x-z);
      [f] = oracleObjectiveModel2(nbObservations,nbArrivals,z,C,D,T,R,regressor,nbRegressors,durations,indexBeta);
      fprintf("k = %i, fold = %.5f, rhs = %.5f, beta_k = %.10f, f = %.5f\n",k,fold, rhs, betak, ...
          f);
      while (bool==1)
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel2(nbObservations,nbArrivals,zAux,C,D,T,R,regressor,nbRegressors,durations,indexBeta);        
        fprintf("\tj = %i, f = %.5f 1/2j = %.5f\n",j, f, 1/(2^j))
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      %upper=min(upper,f);
      betak=bparam/(2^j);
      
      % input("");
      k=k+1;
end
