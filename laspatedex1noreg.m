
function [lambda,obj]=laspatedex1noreg()

cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\LASPATED';

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
C=1;
nbWeeks=20;
nbYears=1;
nbObervationsTotal=nbWeeks*nbYears;
maxObs=nbObervationsTotal;
lambdaTeoricos=zeros(T,R,C);
durations=ones(1,T);
nbVars=1;
sample=zeros(T,R,C,maxObs);
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
nbObservations=zeros(C,R,T);
for t=1:T
    for i=1:R
        for p=1:C
            %lambdaTeoricos(t,g,i,p)=2*rand*duration(t);
            nbObservations(p,i,t)=nbObervationsTotal;
            sample(t,i,p,:)=poissrnd(lambdaTeoricos(t,i,p)*durations(t),1,nbObervationsTotal);
            nbArrivals(p,i,t)=sum(sample(t,i,p,:));
            estimated(t,i,p)=nbArrivals(p,i,t)/(nbObservations(p,i,t)*durations(t));
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
            est(1,t)=est(1,t)+nbArrivals(1,i,t);
            nbtype1=nbtype1+nbObervationsTotal;
        else
            est(2,t)=est(2,t)+nbArrivals(1,i,t);
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
lambda0=epsilon*ones(C,R,T);
weight=1;
sigma=0.5;
betaBar=1;
iterMax=30;
model='noreg';
params.neighbors=neighbors;
params.type=type;
params.distance=distance;
params.alpha=1;
[lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight,sigma,iterMax,epsilon,lambda0,params);
proportion=0.2;
[cputime,weight,lambda]=crossValidation(model,sample,T,R,P,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);


