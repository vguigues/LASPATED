
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
nbWeeks=20;
nbYears=1;
nbObervationsTotal=nbWeeks*nbYears;
maxObs=nbObervationsTotal;
durations=ones(1,T);
nbVars=1;
sample=zeros(T,R,C,nbObervationsTotal);
compteur=1;

for k=1:length(Groups{1,1})
            Theoreticallambdar(Groups{1,1}(k),1)=0.5;
end
for k=1:length(Groups{1,2})
            Theoreticallambdar(Groups{1,2}(k),1)=0.1;
end
for k=1:length(Groups{1,3})
            Theoreticallambdar(Groups{1,3}(k),1)=0.5;
end
for k=1:length(Groups{1,4})
            Theoreticallambdar(Groups{1,4}(k),1)=0.1;
end

for k=1:length(Groups{1,1})
            Theoreticallambdab(Groups{1,1}(k),1)=0.1;
end
for k=1:length(Groups{1,2})
            Theoreticallambdab(Groups{1,2}(k),1)=0.5;
end
for k=1:length(Groups{1,3})
            Theoreticallambdab(Groups{1,3}(k),1)=0.1;
end
for k=1:length(Groups{1,4})
            Theoreticallambdab(Groups{1,4}(k),1)=0.5;
end


for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
end

rmap=[1,2,3,4,5,11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45,56,57,58,59,60,66,67,68,69,70,76,77,78,79,80,86,87,88,89,90,96,97,98,99,100];
bmap=[6,7,8,9,10,16,17,18,19,20,26,27,28,29,30,36,37,38,39,40,46,47,48,49,50,51,52,53,54,55,61,62,63,64,65,71,72,73,74,75,81,82,83,84,85,91,92,93,94,95];

nbArrivals=zeros(C,R,T);
nbObservations=nbObervationsTotal*ones(C,R,T);
for t=1:T
        for p=1:C
            thissample=poissrnd(Theoreticallambdar(t,p)*durations(t),1,nbObervationsTotal);
            for i=1:nbObervationsTotal
                samplesites=1+floor(50*rand(1,thissample(i)));
                for j=1:thissample(i)
                    if (samplesites(j)==51)
                        samplesites(j)=50;
                    end
                    sample(t,rmap(samplesites(j)),p,i)=sample(t,rmap(samplesites(j)),p,i)+1;
                    nbArrivals(p,rmap(samplesites(j)),t)=nbArrivals(p,rmap(samplesites(j)),t)+1;
                end
            end
            thissample=poissrnd(Theoreticallambdab(t,p)*durations(t),1,nbObervationsTotal);
            for i=1:nbObervationsTotal
                samplesites=1+floor(50*rand(1,thissample(i)));
                for j=1:thissample(i)
                    if (samplesites(j)==51)
                        samplesites(j)=50;
                    end
                    sample(t,bmap(samplesites(j)),p,i)=sample(t,bmap(samplesites(j)),p,i)+1;
                    nbArrivals(p,bmap(samplesites(j)),t)=nbArrivals(p,bmap(samplesites(j)),t)+1;
                end
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
%params.alpha=alphas(it);
params.weight=weight*ones(1,nbGroups);   
[lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);








weights=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
alphas=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

err=[];


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


for t=1:T
    for i=1:R
        for p=1:C
            estimated(t,i,p)=nbArrivals(p,i,t)/(nbObservations(p,i,t)*durations(t));
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

est=zeros(2,4);
for t=1:4
    nbtype1=0;
    nbtype2=0;
    for i=1:R
        if (type(i)==1)
            for k=1:7
                est(1,t)=est(1,t)+nbArrivals(1,i,4*(k-1)+t);
                nbtype1=nbtype1+nbObervationsTotal;
            end    
        else
            for k=1:7
                est(2,t)=est(2,t)+nbArrivals(1,i,4*(k-1)+t);
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


for it=1:length(weights)
    it
    params.alpha=alphas(it);
    %params.alpha=0;
    weight=weights(it);
    %weight=220;
    %params.alpha=220;
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight*ones(1,nbGroups),sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambda(p,i,t))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambda(p,i,t))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
    err=[err,mean(differenceL2)]; 
end


plot(weights,err,'k-',LineWidth=1.5);
hold on;
plot(weights,err1,'-.r',LineWidth=1.5);
hold on;
plot(weights,err2,':b',LineWidth=1.5);
hold on;
plot(weights,err3,'--k',LineWidth=1.5);
hold on;
plot(weights,0.028*ones(1,length(weights)),'-r',LineWidth=1.5)

[a,b]=min(err);
hold on;
plot(weights(b),a,'kd',LineWidth=1.5);
[a1,b1]=min(err1);
hold on;
plot(weights(b1),a1,'-.rd',LineWidth=1.5);
hold on;
[a2,b2]=min(err2);
plot(weights(b2),a2,':bd',LineWidth=1.5);
[a3,b3]=min(err3);
hold on;
plot(weights(b3),a3,'--kd',LineWidth=1.5);

legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors';'Cross validation       ']);

differenceL2=[];
for t=1:4
    for p=1:C
           differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-est(1,t))/Theoreticallambdar(t,p)];
           differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-est(2,t))/Theoreticallambdab(t,p)];
    end
end
mean(differenceL2);
est(1,t)
est(2,t)


xT=[];
xr=[];
xemp=[];
yT=[];
yr=[];
yemp=[];

for t=1:T
    xT=[xT,[t-1:0.1:t-1+0.9]];
    xr=[xr,[t-1:0.1:t-1+0.9]];
xemp=[xemp,[t-1:0.1:t-1+0.9]];

yT=[yT,Theoreticallambdab(t,1)*ones(1,10)];
yr=[yr,lambda(1,6,t)*ones(1,10)];
yemp=[yemp,estimated(t,6,1)*ones(1,10)];      
end
plot(xT,yT,'-k');
hold on;
plot(xr,yr,'--b');
hold on;
plot(xemp,yemp,'-.r');
legend(['Theoretical          ';'Regularized estimator';'Empirical estimator  '])

Theoreticallambdar(t,p)
lambda(p,i,t)
estimated(t,i,p)


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

lw=length(weights);
weights=weights(5:lw);
proportion=0.2;
[cputime,weight,lambda]=crossValidation(model,sample,T,R,C,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);

s1=[0:0.05:10];
s2=[0:0.05:10];
ls=length(s1);
z=zeros(ls,ls);
ureds=rand(1,20);
ublues=rand(1,20);

for i=1:ls
    for j=1:ls
        if  (((s1(i)>=0)&&(s1(i)<=5)&&(s2(j)>5)&&(s2(j)<=10))||((s1(i)>5)&&(s1(i)<=10)&&(s2(j)>=0)&&(s2(j)<=5)))
             for k=1:10
                 z(i,j)=z(i,j)+(1/(2^k))*(ublues(2*k-1)*abs(sin(2*pi*k*s1(i)/10))+ublues(2*k)*abs(sin(2*pi*k*s2(j)/10)));
             end
        else
            for k=1:10
                z(i,j)=z(i,j)+(1/(2^k))*((ureds(2*k-1)+1)*abs(sin(2*pi*k*s1(i)/10))+(ureds(2*k)+1)*abs(sin(2*pi*k*s2(j)/10)));
            end
        end
    end
end

surfl(s1,s2,z);
colormap(pink);
shading interp; 

s1=[0:0.05:10];
s2=[0:0.05:10];
ls=length(s1);
z=zeros(ls,ls);
%odd times 

for i=1:ls
    for j=1:ls
        if  (((s1(i)>=0)&&(s1(i)<=5)&&(s2(j)>5)&&(s2(j)<=10))||((s1(i)>5)&&(s1(i)<=10)&&(s2(j)>=0)&&(s2(j)<=5)))
                 z(i,j)=0.1;
        else
                z(i,j)=0.5
        end
    end
end

surfl(s1,s2,z);
colormap(pink);
shading interp; 
zlim([-0.2 0.6])

s1=[0:0.05:10];
s2=[0:0.05:10];
ls=length(s1);
z=zeros(ls,ls);

trueint=zeros(10,10);

for r=1:100
    j=floor((r-1)/10)+1;
    i=mod(r,10)
    if (i==0)
        i=10;
    end
    x=i-0.5;
    y=j-0.5;
    if (isblue(r))
        trueint(i,j)=5*(x+y);
    else
        trueint(i,j)=x+y;
    end
end


for i=1:10
    for j=1:10
    end
end



%even times 

%1 2 3 4 28

%estimated(t,i,1)

xvalue=[0.5:1:9.5];
yvalue=[0.5:1:9.5];
zemp=zeros(10,10);
eventimestep=2;
for i=1:10
    for j=1:10
        zemp(i,j)=estimated(eventimestep,10*(j-1)+i,1);
    end
end


for i=1:ls
    for j=1:ls
        if  (((s1(i)>=0)&&(s1(i)<=5)&&(s2(j)>5)&&(s2(j)<=10))||((s1(i)>5)&&(s1(i)<=10)&&(s2(j)>=0)&&(s2(j)<=5)))
                 z(i,j)=s1(i)+s2(j);
        else
                z(i,j)=5*(s1(i)+s2(j));
        end
    end
end
surfl(s1,s2,z);
colormap(pink);
shading interp; 
zlim([-0.2 0.6])
hold on
surf(xvalue,yvalue,zemp);

sample=zeros(T,R,C,nbObervationsTotal);
for t=1:T
    if (mod(t,2)==0)
        lambdaBarB=75;
        LambdaB=625;
        lambdaBarR=20;
        LambdaR=250;
        for n=1:nbObervationsTotal
            N=poissrnd(LambdaB);
            for k=1:N
                contde=1;
                while contde
                    contd=1;
                    while contd
                        x=10*rand;
                        y=10*rand;
                        if  (((x<=5)&&(y>=5)&&(y<=10))||  ((x>=5)&&(x<=10)&&(y<=5)))
                            contd=0;
                        end
                    end
                    us=lambdaBarB*rand;
                    if (us<=5*(x+y))
                        contde=0;
                        if (x==10)
                            indx=10;
                        else
                            indx=floor(x)+1;
                        end
                        if (y==10)
                            indy=10;
                        else
                            indy=floor(y)+1;
                        end
                        zonenb=(indy-1)*10+indx;
                        sample(t,zonenb,1,n)=sample(t,zonenb,1,n)+1;
                    end
                end
            end
            N=poissrnd(LambdaR);
            for k=1:N
                contde=1;
                while contde
                    contd=1;
                    while contd
                        x=10*rand;
                        y=10*rand;
                        if  (((x<=5) && (y<=5))||((x>=5)&&(x<=10)&&(y>=5)&&(y<=10)))
                            contd=0;
                        end
                    end
                    us=lambdaBarR*rand;
                    if (us<=x+y)
                        contde=0;
                        if (x==10)
                            indx=10;
                        else
                            indx=floor(x)+1;
                        end
                        if (y==10)
                            indy=10;
                        else
                            indy=floor(y)+1;
                        end
                        zonenr=(indy-1)*10+indx;
                        sample(t,zonenr,1,n)=sample(t,zonenr,1,n)+1;
                    end
                end
            end
        end
    else
        lambdaBarB=15;
        LambdaB=125;
        lambdaBarR=100;
        LambdaR=1250;
        for n=1:nbObervationsTotal
            N=poissrnd(LambdaB);
            for k=1:N
                contde=1;
                while contde
                    contd=1;
                    while contd
                        x=10*rand;
                        y=10*rand;
                        if  (((x<=5)&&(y>=5)&&(y<=10))||  ((x>=5)&&(x<=10)&&(y<=5)))
                            contd=0;
                        end
                    end
                    us=lambdaBarB*rand;
                    if (us<=x+y)
                        contde=0;
                        if (x==10)
                            indx=10;
                        else
                            indx=floor(x)+1;
                        end
                        if (y==10)
                            indy=10;
                        else
                            indy=floor(y)+1;
                        end
                        zonenb=(indy-1)*10+indx;
                        sample(t,zonenb,1,n)=sample(t,zonenb,1,n)+1;
                    end
                end
            end

            N=poissrnd(LambdaR);
            for k=1:N
                contde=1;
                while contde
                    contd=1;
                    while contd
                        x=10*rand;
                        y=10*rand;
                        if  (((x<=5) && (y<=5))||((x>=5)&&(x<=10)&&(y>=5)&&(y<=10)))
                            contd=0;
                        end
                    end
                    us=lambdaBarR*rand;
                    if (us<=5*(x+y))
                        contde=0;
                        if (x==10)
                            indx=10;
                        else
                            indx=floor(x)+1;
                        end
                        if (y==10)
                            indy=10;
                        else
                            indy=floor(y)+1;
                        end
                        zonenr=(indy-1)*10+indx;
                        sample(t,zonenr,1,n)=sample(t,zonenr,1,n)+1;
                    end
                end
            end
        end
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


err=[];

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
                    differenceL2=[differenceL2;abs(Theoreticallambda(t,r)-lambda(p,i,t))/Theoreticallambda(t,r)];
                end
            end
    end
    err=[err,mean(differenceL2)]; 
end

plot(weights,err,'k-',LineWidth=1.5);

for t=1:T
    for i=1:R
        for p=1:C
            estimated(t,i,p)=nbArrivals(p,i,t)/(nbObservations(p,i,t)*durations(t));
        end
    end
end

for t=1:T
    for i=1:R
         estimated(t,i,1)
         Integratedlambda(t,i)
         pause
    end
end






differenceL2=[];
for t=1:T
    for i=1:R
        for p=1:C
            differenceL2=[differenceL2;abs(Theoreticallambda(t,r)-estimated(t,i,p))/Theoreticallambda(t,r)];
        end
    end
end
mean(differenceL2)







