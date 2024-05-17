                     
%x=[beta_{c d t}]

%y=[beta_{c d t}]

%min_{y} 0.5 y^T y - x^T y
%beta_{c d t}^T regressor(:,i) >=epsilon   for all c,d,t

%beta_{c d t}=[beta_{c d t 1};beta_{c d t 2};...;beta_{c d t nb_Land_Types+1}]
                                
function [x]=projectionRegressor(regressor,indexBeta,nbRegressors,C,D,T,R,epsilon,sizex,x,lbounds)

model.Q=sparse([0.5*eye(sizex)]);
model.obj=-x';
model.rhs=epsilon*ones(1,C*D*T*R);
model.A=sparse(C*D*T*R,sizex);
nb=1;
ubounds=(10^5)*ones(sizex,1);
for c=1:C
    for d=1:D
        for t=1:T
            ubounds(indexBeta(c,d,t,1))=1;     
            for i=1:R
                for j=1:nbRegressors
                    model.A(nb,indexBeta(c,d,t,j))=regressor(j,i);
                end
                nb=nb+1;
            end
        end
    end
end
model.lb=lbounds;
model.ub=ubounds;
model.sense = '>';
params.outputflag = 0;
results = gurobi(model,params);
x=results.x;

