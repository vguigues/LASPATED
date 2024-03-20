
function [weight,lambda,obj]=laspatedex2noreg()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial discretization for a tes t problem in the plane:
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
nbGroups=4;
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
nbObservations=nbObervationsTotal*ones(C,R,T);
maxObs=nbObervationsTotal;
durations=ones(1,T);
sample=zeros(T,R,C,maxObs);
compteur=1;

for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        %type(compteur)=1;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        %type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
end

u=0;
isblue=zeros(100,1);
for k=1:5 
    u=u+5;
    for i=1:5
        isblue(u+i)=1;
    end
    u=u+5;
end
for k=1:5
    for i=1:5
        isblue(u+i)=1;
    end
    u=u+10;
end

for i=1:R
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
        %Setting the distance to one is equivalent to not considering
        %distances in penalizations. Rather we use the region type filled
        %above.
        distance(i,j)=1;
    end
end

sample=zeros(T,R,C,nbObervationsTotal);
                    
intensities=zeros(T,100);
for r=1:100
    j=floor((r-1)/10)+1;
    i=mod(r,10);
    if (i==0)
        i=10;
    end
    x=i-0.5;
    y=j-0.5;
    for t=1:T
        if (mod(t,2)==0)
            if (isblue(r))
                trueint=5*(x+y);
            else
                trueint=x+y;
            end
            for n=1:nbObervationsTotal
                sample(t,r,1,n)=poissrnd(trueint,1,1);
            end
        else
            if (isblue(r))
                trueint=x+y;
            else
                trueint=5*(x+y);
            end
            for n=1:nbObervationsTotal
                sample(t,r,1,n)=poissrnd(trueint,1,1);
            end
        end
        intensities(t,r)=trueint;
    end
end

nbArrivals=zeros(1,R,T);
for i=1:R
    for t=1:T
        for n=1:nbObervationsTotal
            nbArrivals(1,i,t)=nbArrivals(1,i,t)+sample(t,i,1,n);
        end
    end
end


Theoreticallambda=zeros(T,R);
Integratedlambda=zeros(T,R);
for r=1:100
    j=floor((r-1)/10)+1;
    i=mod(r,10);
    if (i==0)
        i=10;
    end
    for t=1:T
        if (mod(t,2)==0)
            if (isblue(r))       
                Theoreticallambda(t,r)=5*(i+j-1);
                Integratedlambda(t,r)=2.5*(i^2-(i-1)^2)+2.5*(j^2-(j-1)^2);
            else
                Theoreticallambda(t,r)=i+j-1;
                Integratedlambda(t,r)=0.5*(i^2-(i-1)^2)+2.5*(j^2-(j-1)^2);
            end
        else
            if (isblue(r))       
                Theoreticallambda(t,r)=i+j-1;
                Integratedlambda(t,r)=0.5*(i^2-(i-1)^2)+2.5*(j^2-(j-1)^2);
            else
                Theoreticallambda(t,r)=5*(i+j-1);
                Integratedlambda(t,r)=2.5*(i^2-(i-1)^2)+2.5*(j^2-(j-1)^2);
            end
        end
    end
end

est=zeros(T,R);
for t=1:T
    for r=1:R
        est(t,r)=nbArrivals(1,r,t)/nbObervationsTotal;
        est(t,r)
        Theoreticallambda(t,r)
    end
end


epsilon=0.001;
lambda0=epsilon*ones(C,R,T);
sigma=0.5;
betaBar=1;
iterMax=30;
model='noreg';
params.neighbors=neighbors;
params.type=type;
params.distance=distance;
params.alpha=1;

weights=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
alphas=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

err=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err=[err,mean(differenceL2)]; 
end

err1=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    params.alpha=0;
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err1=[err1,mean(differenceL2)]; 
end

Groups=cell(1,2);
nbGroups=2;
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,2)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
for i=1:2
    for j=1:13
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+2];
    end
end

err2=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err2=[err2,mean(differenceL2)]; 
end


err3=[];
for it=1:length(weights)
    it
    %params.alpha=alphas(it);
    params.alpha=0;
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err3=[err3,mean(differenceL2)]; 
end


for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        %type(compteur)=2;
        type(compteur)=1;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        %type(compteur)=2;
        type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
end

Groups=cell(1,4);
nbGroups=4;
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

err4=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    %params.alpha=0;
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err4=[err4,mean(differenceL2)]; 
end



Groups=cell(1,2);
nbGroups=2;
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,2)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
for i=1:2
    for j=1:13
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+2];
    end
end

err5=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    %params.alpha=0;
    weight=weights(it);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
                end
            end
    end
    err5=[err5,mean(differenceL2)]; 
end

% final=length(weights);
% final=4;
% plot(weights(1:final),err(1:final),'k-',LineWidth=1.5);
% hold on;
% plot(weights(1:final),err4(1:final),'-.r',LineWidth=1.5);
% hold on;
% plot(weights(1:final),err1(1:final),':b',LineWidth=1.5);
% hold on;
% plot(weights(1:final),err2(1:final),'--k',LineWidth=1.5);
% hold on;
% plot(weights(1:final),err5(1:final),'-r',LineWidth=1.5);
% hold on;
% plot(weights(1:final),err3(1:final),'-.k',LineWidth=1.5);
% hold on;
% plot(weights(1:final),0.2*ones(1,final),'-b',LineWidth=1.5);

final=length(weights);
%final=6;
plot(weights(1:final),err(1:final),'k-',LineWidth=1.5);
%hold on;
plot(weights(1:final),err4(1:final),'-.r',LineWidth=1.5);
hold on;
plot(weights(1:final),err1(1:final),':b',LineWidth=1.5);
hold on;
%plot(weights(1:final),err2(1:final),'--k',LineWidth=1.5);
%hold on;
plot(weights(1:final),err5(1:final),'-r',LineWidth=1.5);
hold on;
plot(weights(1:final),err3(1:final),'-.k',LineWidth=1.5);
hold on;
plot(weights(1:final),0.043*ones(1,final),'-b',LineWidth=1.5);

%err 4 neighb color type 1 2 N1
%err4 4 neighb no color type 1 1 N2
%err1 4 no neighb
%err2  2 neighb col type 1 2
%err5 2 neighb no color type 1 1 
%err3  2 no neigh

%legend(['4  groups, neighbors 1 ';'4  groups, neighbors 2 ';'4  groups, no neighbors';'2  groups, neighbors 1 ';'2  groups, neighbors 2 ';'2  groups, no neighbors';'Cross validation       ']);
legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors';'Cross validation       ']);

%lw=length(weights);
%weights=weights(5:lw);
proportion=0.2;
[cputime,weight,lambda]=crossValidation(model,sample,T,R,C,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);
differenceL2=[];
for t=1:T
    for i=1:R
        for p=1:C
            differenceL2=[differenceL2;abs(Theoreticallambda(t,i)-lambda(p,i,t))/Theoreticallambda(t,i)];
        end
   end
end
mean(differenceL2)

%n 1 neighb of same color
%n 2 no color


