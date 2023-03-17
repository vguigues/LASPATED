                     
%x=[beta_{c d t}]

%y=[beta_{c d t}]

%min_{y} 0.5 y^T y - x^T y
%beta_{c d t}^T regressor(:,i) >=epsilon   for all c,d,t

%beta_{c d t}=[beta_{c d t 1};beta_{c d t 2};...;beta_{c d t nb_Land_Types+1}]

function [x]=projectionRegressor(regressor,indexBetaA,nbLandTypes,C,D,T,R,epsilon,sizex,x,lbounds,gnb,cnb,nbC,Groups)

model.Q=sparse([0.5*eye(sizex)]);
model.obj=-x';
model.A=sparse(nbC,sizex);
nb=1;
ubounds=(10^5)*ones(sizex,1);

for m=1:length(Groups{1,gnb})
    d=Groups{1,gnb}(m).d;
    t=Groups{1,gnb}(m).t;
    ubounds(indexBetaA(cnb,d,t,1))=1;
    for i=1:R
        for j=1:(1+nbLandTypes)
            model.A(nb,indexBetaA(cnb,d,t,j))=regressor(j,i);
        end
        nb=nb+1;
    end
end
model.rhs=epsilon*ones(1,nbC);
model.lb=lbounds;
model.ub=ubounds;
model.sense = '>';
params.outputflag = 0;
results = gurobi(model,params);
%results = gurobi(model);
x=results.x;

