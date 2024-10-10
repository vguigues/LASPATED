discretization_dir = "Rect10x10";
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance,sample_missing_calls,nb_missing_calls]=read_calls_reg(discretization_dir);


S=10;
C=P;
sigma=0.5;
iterMax=10;
epsilon=0.0001;
x=epsilon*ones(T,G,R,P);
prob=epsilon*ones(P,G,T);
alpha=1;
durations=0.5*ones(G,T);
nbGroups=8;
Groups=cell(1,8);
whichgroup=zeros(T,G);
weights=ones(1,nbGroups);
for i=1:8
    Groups{1,i}=[];
end

for g=1:5
    for t=13:20
        str.g=g;
        str.t=t;
        whichgroup(t,g)=1;
        Groups{1,1}=[Groups{1,1},str];
    end
end


for g=1:5
    for t=21:36
        str.g=g;
        str.t=t;
        whichgroup(t,g)=2;
        Groups{1,2}=[Groups{1,2},str];
    end
end

for g=1:5
    for t=37:44
        str.g=g;
        str.t=t;
        whichgroup(t,g)=3;
        Groups{1,3}=[Groups{1,3},str];
    end
end

for g=1:4
    for t=45:48
        str.g=g;
        str.t=t;
        whichgroup(t,g)=4;
        Groups{1,4}=[Groups{1,4},str];
    end
end

for t=45:48
      str.g=7;
      str.t=t;
      whichgroup(t,7)=4;
      Groups{1,4}=[Groups{1,4},str];
end

for g=1:5
    for t=1:12
        str.g=g;
        str.t=t;
        whichgroup(t,g)=4;
        Groups{1,4}=[Groups{1,4},str];
    end
end

for t=45:48
    str.g=5;
    str.t=t;
    whichgroup(t,5)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=45:48
    str.g=6;
    str.t=t;
    whichgroup(t,6)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=1:12
    str.g=6;
    str.t=t;
    whichgroup(t,6)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=1:12
    str.g=7;
    str.t=t;
    whichgroup(t,7)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=13:20
    str.g=6;
    str.t=t;
    whichgroup(t,6)=6;
    Groups{1,6}=[Groups{1,6},str];
end

for t=13:20
    str.g=7;
    str.t=t;
    whichgroup(t,7)=6;
    Groups{1,6}=[Groups{1,6},str];
end


for t=21:36
    str.g=6;
    str.t=t;
    whichgroup(t,6)=7;
    Groups{1,7}=[Groups{1,7},str];
end

for t=21:36
    str.g=7;
    str.t=t;
    whichgroup(t,7)=7;
    Groups{1,7}=[Groups{1,7},str];
end

for t=37:44
    str.g=6;
    str.t=t;
    whichgroup(t,6)=8;
    Groups{1,8}=[Groups{1,8},str];
end

for t=37:44
    str.g=7;
    str.t=t;
    whichgroup(t,7)=8;
    Groups{1,8}=[Groups{1,8},str];
end
probs=(1/R)*ones(1,R);


printf("nbObservations\n");
for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                aux = nbObservations(t,g,r,p);
                if aux ~= nbObservationsG(g)
                    printf("ERROR: nbObservations(%d %d %d %d) = %d Is different from %d\n",t,g,r,p,aux, nbObservationsG(g));
                    input("");
                end
            end
        end
    end
end


[x,fVal]=projectedGradientMissingLocation(nbObservations,nbCalls,
    nb_missing_calls,neighbors,G,T,R,C,sigma,x,prob,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights);

fVal