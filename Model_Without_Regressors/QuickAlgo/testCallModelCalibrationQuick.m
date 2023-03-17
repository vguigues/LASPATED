
%File generator.cpp in C++

addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_Without_Regressors\QuickAlgo';

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

xmax=100;
ymax=100;
Nx=10;
Ny=10;
R=Nx*Ny;
T=48;
G=7;
P=3;
nbWeeks=52;
nbYears=2;
nbObervationsG=nbWeeks*nbYears*ones(1,G);
lambdaTeoricos=zeros(T,G,R,P);
duration=0.5*ones(1,T);
sampleCalls=zeros(T,G,R,P,nbYears*nbWeeks);

for t=1:T
    for g=1:G
        for i=1:R
            for p=1:P
                lambdaTeoricos(t,g,i,p)=2*rand*duration(t);
                nbObservations(t,g,i,p)=nbWeeks*nbYears;
                sampleCalls(t,g,i,p,:)=poissrnd(lambdaTeoricos(t,g,i,p),1,nbObervationsG(g));
                nbCalls(t,g,i,p)=sum(sampleCalls(t,g,i,p,:));
                estimated(t,g,i,p)=nbCalls(t,g,i,p)/nbObservations(t,g,i,p);
            end
        end
    end
end

for i=1:Nx
    absc(i)=xmax/(2*Nx)+(xmax/Nx)*(i-1);
end
for i=1:Ny
    ord(i)=ymax/(2*Ny)+(ymax/Ny)*(i-1);
end

for i=1:Nx
     for j=1:Ny
         plot(absc(i),ord(j),'rd');
         hold on;
     end
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
    for j=1:length(neighbors{1,i})
        Xj=mod(neighbors{1,i}(j),Nx);
        if (Xj==0)
          Yj=neighbors{1,i}(j)/Nx;
        else
          Yj=1+((neighbors{1,i}(j)-Xj)/Nx);
        end
        if (Xj==0)
            Xj=Nx;
        end
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

epsilon=10^(-2);
sigma=0.5;
betaBar=1;
iterMax=10;
alpha=0.01;
maxDiff=-inf;
meanDiff=0;
%First projected gradient
finalx=zeros(T,G,R,P);
for t=1:T
    for g=1:G
        for p=1:P
            %We compute lambda_{t g p}
            x=epsilon*ones(R,1);
            [lambda1,fVal1]=projectedGradientArmijoBoundaryQuick(nbObservations,nbCalls,neighbors,type,distance,R,sigma,betaBar,x,iterMax,alpha,epsilon,t,g,p);
            for i=1:R
                finalx(t,g,i,p)=lambda1(i);
                difference=abs(lambdaTeoricos(t,g,i,p)-lambda1(i));
                maxDiff=max(maxDiff,difference);
                meanDiff=meanDiff+difference;
            end
        end
    end
end
meanDiff=meanDiff/(T*G*R*P);

epsilon=10^(-2);
finalx2=zeros(T,G,R,P);
meanDiff2=0;

for t=1:T
    for g=1:G
        for p=1:P
            %We compute lambda_{t g p}
            x=epsilon*ones(R,1);
            [lambda1,fVal1]=projectedGradientArmijoFeasibleQuick(nbObservations,nbCalls,neighbors,type,distance,R,sigma,x,iterMax,alpha,epsilon,t,g,p);
            for i=1:R
                finalx2(t,g,i,p)=lambda1(i);
                difference=abs(lambdaTeoricos(t,g,i,p)-lambda1(i));
                meanDiff2=meanDiff2+difference;
            end
        end
    end
end
meanDiff2=meanDiff2/(T*G*R*P);

%Cross validation
proportion=0.2;
x=epsilon*ones(R,1);
[alpha,lambda]=crossValidationQuick(nbObervationsG,sampleCalls,neighbors,type,distance,T,G,R,P,sigma,betaBar,x,iterMax,proportion,epsilon);
