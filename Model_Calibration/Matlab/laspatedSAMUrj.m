

function [lambda,lambdacov]=laspatedSAMUrj()

cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rectangular discretization 10x10
%No regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%discretization_dir='C:\Users\vince\Dropbox\Softwares\LASPATED\Data\Rect10x10_km';
discretization_dir='C:\Users\vince\Dropbox\Softwares\LASPATED\Data\Rect10x10_km';
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls_reg(discretization_dir);

sample_callsA=sample_calls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
C=P;
nbObservations=zeros(P,R,7*T);
nbCalls=zeros(P,R,7*T);
sampleCalls=zeros(7*T,R,P,min(nbObservationsG));
nbObervationsTotal=min(nbObservationsG)*7;

for r=1:R
    for p=1:P
        for g=1:G
            init=(g-1)*T;
            for t=1:T
                nbObservations(p,r,init+t)=nbObservationsA(t,g,r,p);
                nbCalls(p,r,init+t)=nbCallsA(t,g,r,p);
                for j=1:min(nbObservationsG)
                    sampleCalls(init+t,r,p,j)=sample_callsA(t,g,r,p,j);
                end
            end
        end
    end
end

T=7*T;
Groups=cell(1,T);
for t=1:T
    Groups{1,t}=[t];
    whichgroup(t)=t;
end
epsilon=0.001;
lambda0=ones(P,R,T);
sigma=0.5;
betaBar=1;
iterMax=30;
weight=0.1;
durations=0.5*ones(1,T);
x=ones(T,R,P);
nbObservationsTotal=size(sampleCalls,4);
model='noreg';
params.neighbors=neighbors;
params.type=type;
params.distance=distance;
params.alpha=0.1;
nbGroups=length(Groups);
params.weight=weight*ones(1,nbGroups);
[lambda,obj]=laspated(model,nbObservations,nbCalls,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);

weeklambda=zeros(T,1);
for t=1:T
        for r=1:R
            for p=1:P
                weeklambda(t)=weeklambda(t)+lambda(p,r,t);
            end
        end
end
plot(weeklambda,'-k');

heatlambdas=zeros(1,R);
for r=1:R
    for t=1:T
            for p=1:P
                heatlambdas(r)=heatlambdas(r)+x(p,r,t);
            end
    end
end
%dlmwrite('C:\Users\vince\Dropbox\Softwares\LASPATED\Data\heatRect10x10.txt',heatlambdas);
%heatlambdas=dlmread('C:\Users\vince\Dropbox\Softwares\LASPATED\Data\heatRect10x10.txt');

proportion=0.2;
weights=[0,1,10,100,1000,10000];
[cputime,weight,lambda]=crossValidation(model,sampleCalls,T,R,P,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rectangular discretization 10x10
%Regression, no holiday effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

discretization_dir='C:\Users\vince\Dropbox\Softwares\LASPATED\Data\Rect10x10_km';
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls_reg(discretization_dir);
epsilon=10^(-5);
[t1,t2]=size(regressors);
for i=1:t2
    if (sum(regressors(:,i))==0)
        regressors(4,i)=epsilon;
    end
end
C=P;
D=G;
indexBeta=zeros(C,D,T,1+nbLandTypes);
compteur=1;
for c=1:P
    for d=1:G
        for t=1:T
            for j=1:(1+nbLandTypes)
                indexBeta(c,d,t,j)=compteur;       
                compteur=compteur+1;
            end
        end
    end
end
sizex=P*G*T*(1+nbLandTypes);
Groups=cell(1,G*T);
D=G;
compteur=1;
whichgroup=zeros(D,T);
for g=1:G
    for t=1:T
        aux.d=g;
        aux.t=t;
        whichgroup(g,t)=compteur;
        Groups{1,compteur}=[aux];
    end
end

sample_callsA=sample_calls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
C=P;
nbObservations=zeros(P,G,T,R);
nbCalls=zeros(P,G,T,R);
sampleCalls=cell(P,G,T,R);
for r=1:R
    for p=1:P
        for g=1:G
            for t=1:T
                nbObservations(p,g,t,r)=nbObservationsA(t,g,r,p);
                nbCalls(p,g,t,r)=nbCallsA(t,g,r,p);
                for j=1:min(nbObservationsG)
                    sampleCalls{p,g,t,r}=[sampleCalls{p,g,t,r},sample_callsA(t,g,r,p,j)];
                end
            end
        end
    end
end

durations=0.5*ones(1,T);
beta=2*epsilon*ones(sizex,1);
sigma=0.5;
iterMax=15;
nbGroups=length(Groups);
weight=0.1*ones(1,nbGroups);
params.indexBeta=indexBeta;
params.D=D;
params.regressor=regressors;
params.nbRegressors=1+nbLandTypes;
params.sizex=sizex;
%params.nbObs=nbWeeks*7;
%params.isHolidays=isHolidays;
%params.obsbefore=obsbefore;
model='reg';
[lambdacov,obj]=laspated(model,nbObservations,nbCalls,T,R,C,durations,Groups,whichgroup,weight,sigma,iterMax,epsilon,beta,params);

weeklambdaR=zeros(T*G,1);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambdacov(indexBeta(p,g,t,j))*regressors(j,r);
                end
                weeklambdaR(init+t)=weeklambdaR(init+t)+rateest;
            end
        end
    end
end
hold on;
plot(weeklambdaR./(0.5*ones(T*G,1)),'--r');





