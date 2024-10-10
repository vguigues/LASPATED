

function [T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance,sample_missing_calls,nb_missing_calls]=read_calls_reg(discretization_dir)

info_file_path=[discretization_dir '\info.dat'];
calls_file_path=[discretization_dir '\calls.dat'];
neighbors_file_path=[discretization_dir '\neighbors.dat'];
sample_location_file_path=[discretization_dir '\samples.dat'];
missing_calls_file_path=[discretization_dir '\missing_calls.dat'];

%bases_path=[bases_dir '\bases.txt'];
%hospitals_path=[hospitals_dir '\hospitals.txt'];

info_file = fopen(info_file_path);
if (info_file~=-1)
    [line] = strsplit(fgetl(info_file));
    T = str2num(line{1,1});
    G = str2num(line{1,2});
    R = str2num(line{1,3});
    P = str2num(line{1,4});
    nbLandTypes = str2num(line{1,5})-2;
    nbHolidaysYear = str2num(line{1,6});
    line=strsplit(fgetl(info_file));
    nbObservationsG=zeros(1,G);
    for i=1:G
        nbObservationsG(i)=str2num(line{1,i});
    end
    fclose(info_file);
end

calls_file=fopen(calls_file_path);

if (calls_file~=-1)
    sample_calls=zeros(T,G,R,P,max(nbObservationsG));
    nbCalls=zeros(T,G,R,P);
    nbObservations=zeros(T,G,R,P);
    continueBool=1;
    
    while (continueBool)
        str=fgetl(calls_file);
        [line]=strsplit(str);
        if ( (length(str)==3) & (str=='END'))
            continueBool=0;
        else
            t = str2num(line{1,1})+1;
            g = str2num(line{1,2})+1;
            r = str2num(line{1,3})+1;
            p = str2num(line{1,4})+1;
            j = str2num(line{1,5})+1;
            val = str2num(line{1,6});
            h = str2num(line{1,7})+1;
            %if (h==-1)
            sample_calls(t,g,r,p,j) = val;
            nbCalls(t,g,r,p)=nbCalls(t,g,r,p)+sample_calls(t,g,r,p,j);
            %else
            %   sample_calls(t,g,r,p,j) = val;
            %   nbCalls(t,g,r,p)=nbCalls(t,g,r,p)+sample_calls(t,g,r,p,j)
            %end
        end
    end
    fclose(calls_file);
end
%nbObservationsG(8)=30;

estimated=zeros(T,G,R,P);
for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                %if (g<=G)
                    nbObservations(t,g,r,p)=nbObservationsG(g);
                %end
                estimated(t,g,r,p)=nbCalls(t,g,r,p)/nbObservations(t,g,r,p);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

type = zeros(R,1);
regressors = zeros(nbLandTypes+1, R);
neighbors=cell(1,R);
distance = zeros(R,R);

neighbors_file = fopen(neighbors_file_path);
continueBool=1;

while continueBool
    str=fgetl(neighbors_file);
    str=strtrim(str);
    [line]=strsplit (str);
    if (length(str) == 3 & str=='END')
        continueBool=0;
    else
        i=str2num(line{1,1})+1;
        land_type=str2num(line{1,4})+1;
        type(i)=land_type;
        neighbors{1,i}=[];
        for l=1:(1+nbLandTypes)
            regressors(l,i)=str2num(line{1,l+4});
        end
        regressors(1+nbLandTypes,i)=regressors(1+nbLandTypes,i)+str2num(line{1,6+nbLandTypes});
        for l=7+nbLandTypes:2:length(line)
            j=str2num(line{1,l})+1;
            d=str2num(line{1,l+1});
            neighbors{1,i}=[neighbors{1,i},j];
            distance(i,j)=d;
        end
    end
end
fclose(neighbors_file);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
sample_location_file=fopen(sample_location_file_path);
NbLocationSamples=100;
continueBool=1;
while continueBool
    str=fgetl(sample_location_file);
    if ((length(str)==3) & (str=='END'))
        continueBool=0;
    else
        [line]=strsplit (str);
        i=str2num(line{1,1})+1;
        j=str2num(line{1,2})+1;
        sampleLocations(i,j).lat=str2num(line{1,3})+1;
        sampleLocations(i,j).long=str2num(line{1,4})+1;
    end
end

fclose(sample_location_file);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% bases_file=fopen(bases_path);
% continueBool=1;
% i=1;
% while continueBool
%     str=fgetl(bases_file);
%     if ((length(str)==3) & (str=='END'))
%         continueBool=0;
%     else
%         [line]=strsplit (str);
%         bases(i).lat=str2num(line{1,1});
%         bases(i).long=str2num(line{1,2});
%         i=i+1;
%     end
% end
% fclose(bases_file);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% hospitals_file=fopen(hospitals_path);
% continueBool=1;
% i=1;
% while continueBool
%     str=fgetl(hospitals_file);
%     if ((length(str)==3) & (str=='END'))
%         continueBool=0;
%     else
%         [line]=strsplit (str);
%         [line]=strsplit (str);
%         hospitals(i).lat=str2num(line{1,1});
%         hospitals(i).long=str2num(line{1,2});
%         i=i+1;
%     end
% end
% fclose(hospitals_file);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %We assume the missing data is always the location (corresponds
% %to RJ data set).
 missing_calls_file=fopen(missing_calls_file_path);
 if (missing_calls_file~=-1)
     sample_missing_calls=zeros(T,G,P,max(nbObservationsG));
    nb_missing_calls=zeros(T,G,P);
    continueBool=1;
    while (continueBool)
        str=fgetl(missing_calls_file);
        [line]=strsplit(str);
        if ( (length(str)==3) & (str=='END'))
            continueBool=0;
        else
            t = str2num(line{1,1})+1;
            g = str2num(line{1,2})+1;
            p = str2num(line{1,4})+1;
            j = str2num(line{1,5})+1;
            val = str2num(line{1,6});
            h = str2num(line{1,7})+1;
            sample_missing_calls(t,g,p,j)=val;
            nb_missing_calls(t,g,p)=nb_missing_calls(t,g,p)+val;
        end
    end
    fclose(missing_calls_file);
end


