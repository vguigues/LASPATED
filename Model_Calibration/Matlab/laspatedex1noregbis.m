
function [weight,lambda,obj]=laspatedex1noreg()

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
nbWeeks=2;
nbYears=1;
nbObervationsTotal=nbWeeks*nbYears;
maxObs=nbObervationsTotal;
Theoreticallambda=zeros(T,C);
durations=6*ones(1,T);
nbVars=1;
sample=zeros(T,R,C,maxObs);
compteur=1;

for k=1:length(Groups{1,1})
            Theoreticallambdar(Groups{1,1}(k),1)=0.5;
end
for k=1:length(Groups{1,3})
            Theoreticallambdar(Groups{1,3}(k),1)=0.5;
end
for k=1:length(Groups{1,2})
            Theoreticallambdar(Groups{1,2}(k),1)=0.1;
end
for k=1:length(Groups{1,4})
            Theoreticallambdar(Groups{1,4}(k),1)=0.1;
end

for k=1:length(Groups{1,1})
            Theoreticallambdab(Groups{1,1}(k),1)=0.1;
end
for k=1:length(Groups{1,3})
            Theoreticallambdab(Groups{1,3}(k),1)=0.1;
end
for k=1:length(Groups{1,2})
            Theoreticallambdab(Groups{1,2}(k),1)=0.5;
end
for k=1:length(Groups{1,4})
            Theoreticallambdab(Groups{1,4}(k),1)=0.5;
end

for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        for k=1:length(Groups{1,1})
            Theoreticallambda(Groups{1,1}(k),1)=0.5;
        end
        for k=1:length(Groups{1,3})
            Theoreticallambda(Groups{1,3}(k),1)=0.5;
        end
        for k=1:length(Groups{1,2})
            Theoreticallambda(Groups{1,2}(k),1)=0.1;
        end
        for k=1:length(Groups{1,4})
            Theoreticallambda(Groups{1,4}(k),1)=0.1;
        end
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        for k=1:length(Groups{1,1})
            Theoreticallambda(Groups{1,1}(k),1)=0.1;
        end
        for k=1:length(Groups{1,2})
            Theoreticallambda(Groups{1,2}(k),1)=0.5;
        end
        for k=1:length(Groups{1,3})
            Theoreticallambda(Groups{1,3}(k),1)=0.1;
        end
        for k=1:length(Groups{1,4})
            Theoreticallambda(Groups{1,4}(k),1)=0.5;
        end
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        for k=1:length(Groups{1,1})
            Theoreticallambda(Groups{1,1}(k),1)=0.1;
        end
        for k=1:length(Groups{1,2})
            Theoreticallambda(Groups{1,2}(k),1)=0.5;
        end
        for k=1:length(Groups{1,3})
            Theoreticallambda(Groups{1,3}(k),1)=0.1;
        end
        for k=1:length(Groups{1,4})
            Theoreticallambda(Groups{1,4}(k),1)=0.5;
        end
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        for k=1:length(Groups{1,1})
            Theoreticallambda(Groups{1,1}(k),1)=0.5;
        end
        for k=1:length(Groups{1,3})
            Theoreticallambda(Groups{1,3}(k),1)=0.5;
        end
        for k=1:length(Groups{1,2})
            Theoreticallambda(Groups{1,2}(k),1)=0.1;
        end
        for k=1:length(Groups{1,4})
            Theoreticallambda(Groups{1,4}(k),1)=0.1;
        end
        compteur=compteur+1;
    end
end

rmap=[1,2,3,4,5,11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45,56,57,58,59,60,66,67,68,69,70,76,77,78,79,80,86,87,88,89,90,96,97,98,99,100];
bmap=[6,7,8,9,10,16,17,18,19,20,26,27,28,29,30,36,37,38,39,40,46,47,48,49,50,51,52,53,54,55,61,62,63,64,65,71,72,73,74,75,81,82,83,84,85,91,92,93,94,95];

nbObservations=zeros(C,R,T);
for t=1:T
        for p=1:C
            nbObservations(p,i,t)=nbObervationsTotal;
            thissample=poissrnd(Theoreticallambda(t,p)*durations(t)*50,1,nbObervationsTotal);
            for i=1:nbObervationsTotal
                samplesites=1+floor(50*rand(1,thissample(i)));
                for j=1:thissample(i)
                    if (samplesites(j)==51)
                        samplesites(j)=50;
                    end
                    sample(t,rmap(samplesites(j)),p,i)=sample(t,rmap(samplesites(j)),p,i)+1;
                end
            end
            for j=1:nbObervationsTotal
                sample(t,i,p,:)=poissrnd(Theoreticallambda(t,p)*durations(t)*50,1,nbObervationsTotal);
                nbArrivals(p,i,t)=sum(sample(t,i,p,:));
                estimated(t,i,p)=nbArrivals(p,i,t)/(nbObservations(p,i,t)*durations(t));
                nbVars=nbVars+1;
            end
        end
    end
end

%fid=fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\sample6.txt','w');
%for t=1:T
%    for i=1:R
%        for p=1:C
%            for k=1:nbObervationsTotal
%            fprintf(fid,'%d %d %d %d %d\n',[t i p k sample(t,i,p,k)]);
%            end
%        end
%    end
%end
%fclose(fid);

est=zeros(2,2);
for t=1:2
    nbtype1=0;
    nbtype2=0;
    for i=1:R
        if (type(i)==1)
            for k=1:28
                est(1,t)=est(1,t)+nbArrivals(1,i,2*k-1);
                nbtype1=nbtype1+nbObervationsTotal;
            end    
        else
            for k=1:28
                est(2,t)=est(2,t)+nbArrivals(1,i,2*k);
                nbtype2=nbtype2+nbObervationsTotal;
            end
        end
    end
    est(1,t)=est(1,t)/(nbtype1*durations(t));
    est(2,t)=est(2,t)/(nbtype2*durations(t));
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
        %Setting the distance to one is equivalent to not considering
        %distances in penalizations. Rather we use the region type filled
        %above.
        distance(i,j)=1;
    end
end

epsilon=0.001;
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

weights=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
alphas=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

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
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,i,p)-lambda(p,i,t))/Theoreticallambda(t,i,p)];
                end
            end
    end
    err=[err,mean(differenceL2)]; 
end
figure()
plot(err);

differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambda(t,i,p)-est(1,t))/Theoreticallambda(t,i,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambda(t,i,p)-est(2,t))/Theoreticallambda(t,i,p)];
                    end
                end
            end
    end
mean(differenceL2)

proportion=0.2;
[cputime,weight,lambda]=crossValidation(model,sample,T,R,C,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);


