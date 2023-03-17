
%File generator.cpp is the C++ version of the code in this folder

cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\LASPATED';

%Test with real data
discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Rect10x10_km';
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls(discretization_dir);

sample_callsA=sample_calls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
C=P;
nbObservations=zeros(7*T,R,P);
nbCalls=zeros(7*T,R,P);
sampleCalls=zeros(7*T,R,P,104);
%sampleCalls=zeros(7*T,R,P,max(nbObservationsG));
for r=1:R
    for p=1:P
        for g=1:G
            init=(g-1)*T;
            for t=1:T
                nbObservations(init+t,r,p)=nbObservationsA(t,g,r,p);
                %size(nbObservations)
                %init+t
                nbCalls(init+t,r,p)=nbCallsA(t,g,r,p);
                for j=1:104
                    sampleCalls(init+t,r,p,j)=sample_callsA(t,g,r,p,j);
                    %sampleCalls(init+t,r,p,j)=sample_callsA{p,g,t,r}(j);
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
alpha=0;
epsilon=0.001;
x=ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=30;
durations=6*ones(1,T);
weight=0;

%Call model calibration without regressors

%[xnr,fVal]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight)
[x,fVal]=  projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight)
weeklambda=zeros(T,1);
for t=1:T
        for r=1:R
            for p=1:P
                weeklambda(t)=weeklambda(t)+x(t,r,p);
                %weeklambda(t)=weeklambda(t)+estimated(t,g,r,p);
            end
        end
end
plot(weeklambda);

weeklambdaprioritites=zeros(T,P);
for t=1:T
    for r=1:R
        for p=1:P
            weeklambdaprioritites(t,p)=weeklambdaprioritites(t,p)+x(t,r,p);  
        end
    end
end
plot(weeklambdaprioritites(:,1));
hold on;
plot(weeklambdaprioritites(:,2));
hold on;
plot(weeklambdaprioritites(:,3));
legend(['Priority 1';'Priority 2';'Priority 3']);
xnr=x;
x=ones(T,R,P);
[cputime,alpha,weight,lambda]=crossValidation(nbObervationsTotal,sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup);

%Call model calibration with regressors

[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,
    nbCalls,nbObservations,estimated,type,regressors,
    neighbors,distance]=read_calls(discretization_dir);

T=8*T;
%sample_calls(T,G+1,R,P,max(nbObservationsG))
sample_calls(T,G+1,R,P,104)
nbCalls(T,G+1,R,P);
nbObservations(T,G+1,R,P);
type(i)
regressors = zeros(land,i);
neighbors{1,i}
distance(i,j)

[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbCalls,x,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup,lbounds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial discretization for test problem in the plane:
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
% 1,2  2,2  *  *  *  *  *  *  *  *
% 1,1  2,1 3,1 *  *  *  *  *  *  *
% (Xi,Yi) below for region numbered i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax=10;
ymax=10;
Nx=10;
Ny=10;
R=Nx*Ny;
T=4*7;
Groups=cell(1,4);
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,4)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
Groups{1,3}=[3];
Groups{1,4}=[4];
for i=1:4
    for j=1:6
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+4];
    end
end
P=1;
nbWeeks=20;
nbYears=1;
nbObervationsTotal=nbWeeks*nbYears;
maxObs=nbObervationsTotal;
lambdaTeoricos=zeros(T,R,P);
durations=ones(1,T);
nbVars=1;
sampleCalls=zeros(T,R,P,maxObs);
GroupAreas=cell(1,4);
compteur=1;
for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.1;
        end
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.5;
        end
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.5;
        end
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.1;
        end
        compteur=compteur+1;
    end
end

% r=1;
% for i=1:length(Groups{1,1})
%     lambdaTeoricos(Groups{1,1}(i),r,1)=10;
%     lambdaTeoricos(Groups{1,1}(i),4,1)=10;
% end
% for i=1:length(Groups{1,3})
%     lambdaTeoricos(Groups{1,3}(i),r,1)=10;
%     lambdaTeoricos(Groups{1,3}(i),4,1)=10;
% end
% for i=1:length(Groups{1,2})
%     lambdaTeoricos(Groups{1,2}(i),r,1)=2;
%     lambdaTeoricos(Groups{1,2}(i),4,1)=2;
% end
% for i=1:length(Groups{1,4})
%     lambdaTeoricos(Groups{1,4}(i),r,1)=2;
%     lambdaTeoricos(Groups{1,4}(i),4,1)=2;
% end
% 
% r=3;
% for i=1:length(Groups{1,1})
%     lambdaTeoricos(Groups{1,1}(i),r,1)=2;
%     lambdaTeoricos(Groups{1,1}(i),2,1)=2;
% end
% for i=1:length(Groups{1,2})
%     lambdaTeoricos(Groups{1,2}(i),r,1)=10;
%     lambdaTeoricos(Groups{1,2}(i),2,1)=10;
% end
% for i=1:length(Groups{1,3})
%     lambdaTeoricos(Groups{1,3}(i),r,1)=2;
%     lambdaTeoricos(Groups{1,3}(i),2,1)=2;
% end
% for i=1:length(Groups{1,4})
%     lambdaTeoricos(Groups{1,4}(i),r,1)=10;
%     lambdaTeoricos(Groups{1,4}(i),2,1)=10;
% end
% 
for t=1:T
    for i=1:R
        for p=1:P
            %lambdaTeoricos(t,g,i,p)=2*rand*duration(t);
            nbObservations(t,i,p)=nbObervationsTotal;
            sampleCalls(t,i,p,:)=poissrnd(lambdaTeoricos(t,i,p)*durations(t),1,nbObervationsTotal);
            nbCalls(t,i,p)=sum(sampleCalls(t,i,p,:));
            estimated(t,i,p)=nbCalls(t,i,p)/(nbObservations(t,i,p)*durations(t));
            nbVars=nbVars+1;
        end
    end
end

est=zeros(2,4);
for t=1:4
    nbtype1=0;
    nbtype2=0;
    for i=1:R
        if (type(i)==1)
            est(1,t)=est(1,t)+nbCalls(t,i,1);
            nbtype1=nbtype1+nbObervationsTotal;
        else
            est(2,t)=est(2,t)+nbCalls(t,i,1);
            nbtype2=nbtype2+nbObervationsTotal;
        end
    end
    est(1,t)=est(1,t)/nbtype1;
    est(2,t)=est(2,t)/nbtype2;
end


for i=1:Nx
    absc(i)=xmax/(2*Nx)+(xmax/Nx)*(i-1);
end
for i=1:Ny
    ord(i)=ymax/(2*Ny)+(ymax/Ny)*(i-1);
end

for i=1:Nx
      for j=1:Ny
          if (type((j-1)*Nx+i)==1)
             plot(absc(i),ord(j),'rs');
             hold on;
          else
             plot(absc(i),ord(j),'bo');
             hold on;
          end
      end
end

for i=1:R
    %type(i)=1+floor(4*rand);  
    Xi=mod(i,Nx);
    if (Xi==0)
       Yi=i/Nx;
    else
       Yi=1+(i-Xi)/Nx;
    end
    if (Xi==0)
        Xi=Nx;
    end
    neighbors{1,i}=[];
    %Possible neighbors are (Xi,Yi-1),(Xi,Yi+1),
    %(Xi-1,Yi),(Xi+1,Yi),
    %(Xi-1,Yi+1),(Xi-1,Yi-1),
    %(Xi+1,Yi+1), and (Xi+1,Yi-1),
    if (Yi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi];     
    end
    if (Yi+1<=Ny)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi];     
    end
    if (Xi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi-1];     
    end
    if (Xi+1<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi+1];     
    end
    if ((Yi+1<=Ny)&&(Xi-1>=1))
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi-1];     
    end
    if ((Yi-1>=1)&&(Xi-1)>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi-1];     
    end
    if ((Yi+1<=Ny)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi+1];     
    end
    if ((Yi-1>=1)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi+1];     
    end
    for j=1:R
        Xj=mod(j,Nx);
        if (Xj==0)
          Yj=j/Nx;
        else
          Yj=1+((j-Xj)/Nx);
        end
        if (Xj==0)
            Xj=Nx;
        end
        %distance(i,j)=norm([absc(Xi)-absc(Xj);ord(Yi)-ord(Yj)]);
        distance(i,j)=1;
    end
end

%Test of projectedGradientArmijoBoundary
epsilon=0.001;
%x=epsilon*ones(sizex,1);
x=epsilon*ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=100;
%First projected gradient

proportion=0.2;
[cputime,alpha,weight,lambda]=crossValidation(nbObervationsTotal,sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup);

weights=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
alphas=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

err=[];

for it=1:length(weights)
it
    weight=weights(it);
alpha=alphas(it);

[lambda1,fVal1]=projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,1,epsilon/2,durations,Groups,whichgroup,1);
%[lambda1,fVal1]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon/2,durations,Groups,whichgroup,weight);
differenceL=[];
differenceL2=[];
for t=1:T
        for i=1:R
            for p=1:P
                differenceL=[differenceL;abs(lambdaTeoricos(t,i,p)-estimated(t,i,p))/lambdaTeoricos(t,i,p)];
                differenceL2=[differenceL2;abs(lambdaTeoricos(t,i,p)-lambda1(t,i,p))/lambdaTeoricos(t,i,p)];
            end
        end
end
err=[err,mean(differenceL2)];
end
% alpha=1;
% 
% [lambda12,fVal12]=projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,alpha,epsilon/2,durations,Groups,whichgroup,weight);
% %[lambda1,fVal1]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon/2,durations,Groups,whichgroup,weight);
% differenceL22=[];
% for t=1:T
%         for i=1:R
%             for p=1:P
%                 differenceL22=[differenceL22;abs(lambdaTeoricos(t,i,p)-lambda12(t,i,p))/lambdaTeoricos(t,i,p)];
%             end
%         end
% end

nb=20;
errTotal=0.387
err20=[0.3866 0.3856 0.3816 0.3856 0.377 0.3205 0.2432    0.3388    0.3997    0.4090    0.0930    0.0815    0.0733];
err20=[err20,0.0671    0.0623    0.0591    0.0556    0.0530    0.0508    0.0490    0.0479    0.0477    0.0467];
err20=[err20,0.0553    0.0579    0.0680    0.0618    0.0553    0.0499    0.0435    0.0412    0.0380    0.0439];
err20=[err20,0.0452    0.0433    0.0487    0.0591    0.0672    0.0726    0.0950    0.1054    0.1148    0.0975];
err20=[err20,0.1184    0.1132    0.1137    0.1330];

err2 empirique  1.20

err100

err1000

weights=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
plot(weights,err20,'k-');
hold on;
plot(weights,err20,'r-.');
hold on;
plot(weights,err100,'c:');
hold on;
plot(weights,err1000,'b--');
hold on;
legend(['2    periods';'20   periods';'100  periods';'1000 periods'])
xTick = [0,1,5,10,100,200,300,400];
set(gca,'xtick',xTick);

mean(differenceL2)
plot(differenceL);
mean(differenceL)
hold on;
plot(differenceL2);
hold on;
%plot(differenceL22);
legend(['Empirical estimator        ';'Regularized estimator      '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_Without_Regressors';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial discretization for test problem in the plane:
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
% 1,2  2,2  *  *  *  *  *  *  *  *
% 1,1  2,1 3,1 *  *  *  *  *  *  *
% (Xi,Yi) below for region numbered i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax=10;
ymax=10;
Nx=2;
Ny=2;
R=Nx*Ny;
T=4;
Groups=cell(1,4);
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,4)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
Groups{1,3}=[3];
Groups{1,4}=[4];
P=1;
nbWeeks=52;
nbYears=40;
nbObervationsTotal=nbWeeks*nbYears*7;
maxObs=nbObervationsTotal;
lambdaTeoricos=zeros(T,R,P);
durations=6*ones(1,T);
nbVars=1;
sampleCalls=zeros(T,R,P,maxObs);
r=1;
for i=1:length(Groups{1,1})
    lambdaTeoricos(Groups{1,1}(i),r,1)=5;
    lambdaTeoricos(Groups{1,1}(i),4,1)=5;
end
for i=1:length(Groups{1,3})
    lambdaTeoricos(Groups{1,3}(i),r,1)=5;
    lambdaTeoricos(Groups{1,3}(i),4,1)=5;
end
for i=1:length(Groups{1,2})
    lambdaTeoricos(Groups{1,2}(i),r,1)=2;
    lambdaTeoricos(Groups{1,2}(i),4,1)=2;
end
for i=1:length(Groups{1,4})
    lambdaTeoricos(Groups{1,4}(i),r,1)=2;
    lambdaTeoricos(Groups{1,4}(i),4,1)=2;
end


r=3;
for i=1:length(Groups{1,1})
    lambdaTeoricos(Groups{1,1}(i),r,1)=2;
    lambdaTeoricos(Groups{1,1}(i),2,1)=2;
end
for i=1:length(Groups{1,2})
    lambdaTeoricos(Groups{1,2}(i),r,1)=5;
    lambdaTeoricos(Groups{1,2}(i),2,1)=5;
end
for i=1:length(Groups{1,3})
    lambdaTeoricos(Groups{1,3}(i),r,1)=2;
    lambdaTeoricos(Groups{1,3}(i),2,1)=2;
end
for i=1:length(Groups{1,4})
    lambdaTeoricos(Groups{1,4}(i),r,1)=5;
    lambdaTeoricos(Groups{1,4}(i),2,1)=5;
end

for t=1:T
    for i=1:R
        for p=1:P
            indexLambda(t,i,p)=nbVars;
            %lambdaTeoricos(t,g,i,p)=2*rand*duration(t);
            nbObservations(t,i,p)=nbObervationsTotal;
            sampleCalls(t,i,p,:)=poissrnd(lambdaTeoricos(t,i,p)*durations(t),1,nbObervationsTotal);
            nbCalls(t,i,p)=sum(sampleCalls(t,i,p,:));
            estimated(t,i,p)=nbCalls(t,i,p)/(nbObservations(t,i,p)*durations(t));
            nbVars=nbVars+1;
        end
    end
end

for i=1:Nx
    absc(i)=xmax/(2*Nx)+(xmax/Nx)*(i-1);
end
for i=1:Ny
    ord(i)=ymax/(2*Ny)+(ymax/Ny)*(i-1);
end

for i=1:R
    %type(i)=1+floor(4*rand);  
    type(i)=1;
    Xi=mod(i,Nx);
    if (Xi==0)
       Yi=i/Nx;
    else
       Yi=1+(i-Xi)/Nx;
    end
    if (Xi==0)
        Xi=Nx;
    end
    neighbors{1,i}=[];
    %Possible neighbors are (Xi,Yi-1),(Xi,Yi+1),
    %(Xi-1,Yi),(Xi+1,Yi),
    %(Xi-1,Yi+1),(Xi-1,Yi-1),
    %(Xi+1,Yi+1), and (Xi+1,Yi-1),
    if (Yi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi];     
    end
    if (Yi+1<=Ny)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi];     
    end
    if (Xi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi-1];     
    end
    if (Xi+1<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi+1];     
    end
    if ((Yi+1<=Ny)&&(Xi-1>=1))
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi-1];     
    end
    if ((Yi-1>=1)&&(Xi-1)>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi-1];     
    end
    if ((Yi+1<=Ny)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi+1];     
    end
    if ((Yi-1>=1)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi+1];     
    end
    for j=1:R
        Xj=mod(j,Nx);
        if (Xj==0)
          Yj=j/Nx;
        else
          Yj=1+((j-Xj)/Nx);
        end
        if (Xj==0)
            Xj=Nx;
        end
        distance(i,j)=norm([absc(Xi)-absc(Xj);ord(Yi)-ord(Yj)]);
    end
end

%Test of projectedGradientArmijoBoundary
epsilon=0.01;
sizex=T*R*P;
%x=epsilon*ones(sizex,1);
x=epsilon*ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=1000;
alpha=0;
%First projected gradient
w=0;
[lambda1,fVal1]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,
type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight);
differenceL3=[];
differenceL4=[];
for t=1:T
        for i=1:R
            for p=1:P
                differenceL3=[differenceL3;abs(lambdaTeoricos(t,i,p)-estimated(t,i,p))];
                differenceL4=[differenceL4;abs(lambdaTeoricos(t,i,p)-lambda01(t,i,p))];
            end
        end
end
hold on;
mean(differenceL3);
plot(differenceL3);
mean(differenceL4);
hold on;
plot(differenceL4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax=10;
ymax=10;
Nx=2;
Ny=2;
R=Nx*Ny;
T=4;
Groups=cell(1,4);
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,4)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
Groups{1,3}=[3];
Groups{1,4}=[4];
P=1;
nbWeeks=52;
nbYears=2;
nbObervationsTotal=nbWeeks*nbYears*7;
maxObs=nbObervationsTotal;
lambdaTeoricos=zeros(T,R,P);
durations=6*ones(1,T);
nbVars=1;
sampleCalls=zeros(T,R,P,maxObs);
r=1;
for i=1:length(Groups{1,1})
    lambdaTeoricos(Groups{1,1}(i),r,1)=5;
    lambdaTeoricos(Groups{1,1}(i),4,1)=5;
end
for i=1:length(Groups{1,3})
    lambdaTeoricos(Groups{1,3}(i),r,1)=5;
    lambdaTeoricos(Groups{1,3}(i),4,1)=5;
end
for i=1:length(Groups{1,2})
    lambdaTeoricos(Groups{1,2}(i),r,1)=2;
    lambdaTeoricos(Groups{1,2}(i),4,1)=2;
end
for i=1:length(Groups{1,4})
    lambdaTeoricos(Groups{1,4}(i),r,1)=2;
    lambdaTeoricos(Groups{1,4}(i),4,1)=2;
end


r=3;
for i=1:length(Groups{1,1})
    lambdaTeoricos(Groups{1,1}(i),r,1)=2;
    lambdaTeoricos(Groups{1,1}(i),2,1)=2;
end
for i=1:length(Groups{1,2})
    lambdaTeoricos(Groups{1,2}(i),r,1)=5;
    lambdaTeoricos(Groups{1,2}(i),2,1)=5;
end
for i=1:length(Groups{1,3})
    lambdaTeoricos(Groups{1,3}(i),r,1)=2;
    lambdaTeoricos(Groups{1,3}(i),2,1)=2;
end
for i=1:length(Groups{1,4})
    lambdaTeoricos(Groups{1,4}(i),r,1)=5;
    lambdaTeoricos(Groups{1,4}(i),2,1)=5;
end

for t=1:T
    for i=1:R
        for p=1:P
            indexLambda(t,i,p)=nbVars;
            %lambdaTeoricos(t,g,i,p)=2*rand*duration(t);
            nbObservations(t,i,p)=nbObervationsTotal;
            sampleCalls(t,i,p,:)=poissrnd(lambdaTeoricos(t,i,p)*durations(t),1,nbObervationsTotal);
            nbCalls(t,i,p)=sum(sampleCalls(t,i,p,:));
            estimated(t,i,p)=nbCalls(t,i,p)/(nbObservations(t,i,p)*durations(t));
            nbVars=nbVars+1;
        end
    end
end

for i=1:Nx
    absc(i)=xmax/(2*Nx)+(xmax/Nx)*(i-1);
end
for i=1:Ny
    ord(i)=ymax/(2*Ny)+(ymax/Ny)*(i-1);
end

for i=1:R
    %type(i)=1+floor(4*rand);  
    type(i)=1;
    Xi=mod(i,Nx);
    if (Xi==0)
       Yi=i/Nx;
    else
       Yi=1+(i-Xi)/Nx;
    end
    if (Xi==0)
        Xi=Nx;
    end
    neighbors{1,1}=[4];
    neighbors{1,4}=[1];
    neighbors{1,2}=[3];
    neighbors{1,3}=[2];
    for j=1:R
        Xj=mod(j,Nx);
        if (Xj==0)
          Yj=j/Nx;
        else
          Yj=1+((j-Xj)/Nx);
        end
        if (Xj==0)
            Xj=Nx;
        end
        %distance(i,j)=norm([absc(Xi)-absc(Xj);ord(Yi)-ord(Yj)]);
        distance(i,j)=1;
    end
end

%Test of projectedGradientArmijoBoundary
epsilon=0.01;
sizex=T*R*P;
%x=epsilon*ones(sizex,1);
x=epsilon*ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=1000;
alpha=100;
%First projected gradient
w=0;
[lambda1,fVal1]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight);
differenceL3=[];
differenceL4=[];
for t=1:T
        for i=1:R
            for p=1:P
                differenceL3=[differenceL3;abs(lambdaTeoricos(t,i,p)-estimated(t,i,p))];
                differenceL4=[differenceL4;abs(lambdaTeoricos(t,i,p)-lambda1(t,i,p))];
            end
        end
end
hold on;
mean(differenceL3);
plot(differenceL3);
mean(differenceL4);
hold on;
plot(differenceL4);
 legend(['Empirical ';'Calibrated'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Test of projectedGradientArmijoFeasible
%Second projected gradient
alpha=1;
iterMax=200;
x=epsilon*ones(T,G,R,P);
[lambda01,fVal2]=projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight);

alpha=1000;
x=epsilon*ones(T,G,R,P);
[lambda1000,fVal2]=projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight);

meanp=zeros(1,4);
for i=1:4
g=1;
p=1;
for t=1:T
    lambdaTeoricos(t,g,i,p)-lambda1000(t,g,i,p)
    pause
    meanp(i)=meanp(i)+100*abs(lambdaTeoricos(t,g,i,p)-lambda1000(t,g,i,p))/lambdaTeoricos(t,g,i,p);
end
meanp(i)=meanp(i)/4;
end
%Customizing text in xAxis: content and position.
%Customize x and y ranges.



 type(compteur)=1;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.1;
        end
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        for k=1:length(Groups{1,1})
            lambdaTeoricos(Groups{1,1}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,2})
            lambdaTeoricos(Groups{1,2}(k),compteur,1)=0.5;
        end
        for k=1:length(Groups{1,3})
            lambdaTeoricos(Groups{1,3}(k),compteur,1)=0.1;
        end
        for k=1:length(Groups{1,4})
            lambdaTeoricos(Groups{1,4}(k),compteur,1)=0.5;
        end
       
r=6;
p=1;
x=[];
init=0;
for i=1:28
    x=[x,[init:0.1:init+1]];
    init=init+1;
end
yth=[];
for t=1:28
    yth=[yth,lambdaTeoricos(t,r,p)*ones(1,11)];
end
yreg=[];
for t=1:28
    yreg=[yreg,lambda1(t,r,p)*ones(1,11)];
end
yest=[];
for t=1:28
    yest=[yest,estimated(t,r,p)*ones(1,11)];
end
plot(x,yth,'k-');
hold on;
plot(x,yreg,'r--');
hold on;
plot(x,yest,'b-.');
hold on;
ylim([-0.1 0.9])
xTick =1:28;
set(gca,'xtick',xTick);
%yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
%xTickLabel = {{'first';'label'},'second',{'third';'label'},'long fourth tick label ','fifth'};
xTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28'};
for k = 1:length(xTick)
    if (k==1)
        xpos=xTick(1)/2;
    else
        xpos=xTick(k-1)+0.5*(xTick(k)-xTick(k-1));
    end
    text(xpos,yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center')
end

y01=[lambda1(1,r,p)*ones(1,11),lambda1(2,r,p)*ones(1,11),lambda1(3,r,p)*ones(1,11),lambda1(4,r,p)*ones(1,11)];
hold on;
plot(x,y01,'b-.');

y1000=[lambda1000(1,g,r,p)*ones(1,11),lambda1000(2,g,r,p)*ones(1,11),lambda1000(3,g,r,p)*ones(1,11),lambda1000(4,g,r,p)*ones(1,11)];
hold on;
plot(x,y1000,'r--');


differenceL3=[];
for t=1:T
    for g=1:G
        for i=1:R
            for p=1:P
                differenceL3=[differenceL3;abs(lambdaTeoricos(t,g,i,p)-lambda2(t,g,i,p))];
            end
        end
    end
end
mean(differenceL3)
plot(differenceL);

%Cross validation
proportion=0.2;
iterMax=30;
epsilon=0.01;
sizex=T*G*R*P;
%x=epsilon*ones(sizex,1);
x=epsilon*ones(T,G,R,P);
weights=[0,1];
[cputime,alpha,lambda]= crossValidationNoReg(sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights);

Calls=importdata('C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_Without_Regressors\Rect10x10_km\calls_observations.txt');
T=48;
G=7;
R=76;
P=3;
Nb=size(Calls,1);
nbObservations=zeros(T,G,R,P);
nbCalls=zeros(T,G,R,P);
for i=1:Nb
    nbObservations(Calls(i,1)+1,Calls(i,2)+1,Calls(i,3)+1,Calls(i,4)+1)=Calls(i,5);
    nbCalls(Calls(i,1)+1,Calls(i,2)+1,Calls(i,3)+1,Calls(i,4)+1)=Calls(i,6);
end
type=zeros(76,1);
neighbors=cell(1,76);
regressor=zeros(5,76);
distance=zeros(76,76);
for i=1:76
    i
    neighbors{1,i}=[];
    type(i)=Neigh{i,1}(4)+1;
    regressor(:,i)=Neigh{i,1}(5:9)';
    j=10;
    while (j+1<=length(Neigh{i,1}))
        neighbors{1,i}=[neighbors{1,i},floor(Neigh{i,1}(j))+1];
        distance(i,floor(Neigh{i,1}(j))+1)=Neigh{i,1}(j+1);
        j=j+1;
    end
end

nbLandTypes=5;


Lagg=zeros(T,1);
for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                Lagg(t)=Lagg(t)+lambda2(t,g,r,p);
            end
        end
    end
end

Lagg2=zeros(T,G);
for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                Lagg2(t,g)=Lagg2(t,g)+lambda2(t,g,r,p);
            end
        end
    end
end

