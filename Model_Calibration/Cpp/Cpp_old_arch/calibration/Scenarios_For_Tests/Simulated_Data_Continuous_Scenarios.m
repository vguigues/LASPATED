
%Simulate calls in RJ region

for nbScen=1:nbScenarios
    times=sort(t0*0.5+2*rand(1,40));
    for i=1:40
        r=1+floor(76*rand);
        p=1+floor(3*rand);
        calls{1,nbScen}(i).time=times(i);
        calls{1,nbScen}(i).region=r;
        calls{1,nbScen}(i).priority=p;
        calls{1,nbScen}(i).day=g;
        calls{1,nbScen}(i).time_on_scene=0.1;
        latCall=sampleLocations(r,floor(rand*100)+1).lat;
        longCall=sampleLocations(r,floor(rand*100)+1).long;
        calls{1,nbScen}(i).lat=latCall;
        calls{1,nbScen}(i).long=longCall;
        calls{1,nbScen}(i).time_at_hospital=0.2;
        calls{1,nbScen}(i).timeCleaningBase=1;
        calls{1,nbScen}(i).cleaning_needed=0;
        calls{1,nbScen}(i).hosp_needed=1;
        indexhosp=1;
        mindist=mydistance(latCall,longCall,hospitals(1).lat,hospitals(1).long,'sphere',6371);
        for indexh=2:nb_hospitals
            dist=mydistance(latCall,longCall,hospitals(indexh).lat,hospitals(indexh).long,'sphere',6371);
            if (dist<mindist)
                mindist=dist;
                indexhosp=indexh;
            end
        end
        %calls{1,nbScen}(i).ih=indexhosp;
        calls{1,nbScen}(i).ih=1+floor(9*rand);
    end
end
        
fileEstimation=fopen('C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Scenarios_For_Tests\Simulated_Data_Continuous_Scenarios\Random_Hospital\simualtedRealLatLong40.txt','w');
for nbScen=1:nbScenarios
    fprintf(fileEstimation,'%f\n',40);
    for i=1:40
        fprintf(fileEstimation,'%f %f %f %f %f %f %f %f %f %f %f %f\n',calls{1,nbScen}(i).time,calls{1,nbScen}(i).region,calls{1,nbScen}(i).priority,calls{1,nbScen}(i).day,calls{1,nbScen}(i).time_on_scene,calls{1,nbScen}(i).lat,calls{1,nbScen}(i).long,calls{1,nbScen}(i).time_at_hospital,calls{1,nbScen}(i).timeCleaningBase,calls{1,nbScen}(i).cleaning_needed,calls{1,nbScen}(i).hosp_needed,calls{1,nbScen}(i).ih);
    end
end
fclose(fileEstimation);

%Simulate queues of calls in the RJ area.

times=[1;4;7];
for nbScen=1:nbScenarios
    i=1;
    for t=1:3
        for ind=1:10
            r=1+floor(76*rand);
            p=1+floor(3*rand);
            calls{1,nbScen}(i).time=times(t);
            calls{1,nbScen}(i).region=r;
            calls{1,nbScen}(i).priority=p;
            calls{1,nbScen}(i).day=5;
            calls{1,nbScen}(i).time_on_scene=0.1;
            latCall=sampleLocations(r,floor(rand*100)+1).lat;
            longCall=sampleLocations(r,floor(rand*100)+1).long;
            calls{1,nbScen}(i).lat=latCall;
            calls{1,nbScen}(i).long=longCall;
            calls{1,nbScen}(i).time_at_hospital=0.2;
            calls{1,nbScen}(i).timeCleaningBase=1;
            calls{1,nbScen}(i).cleaning_needed=0;
            calls{1,nbScen}(i).hosp_needed=1;
            indexhosp=1;
            mindist=mydistance(latCall,longCall,hospitals(1).lat,hospitals(1).long,'sphere',6371);
            for indexh=2:nb_hospitals
                dist=mydistance(latCall,longCall,hospitals(indexh).lat,hospitals(indexh).long,'sphere',6371);
                if (dist<mindist)
                    mindist=dist;
                    indexhosp=indexh;
                end
            end
            calls{1,nbScen}(i).ih=indexhosp;
            %calls{1,nbScen}(i).ih=1+floor(9*rand);
            i=i+1;
        end
    end
end
        
fileEstimation=fopen('C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Scenarios_For_Tests\Queues_Of_Calls\simualtedQueues.txt','w');
for nbScen=1:nbScenarios
    fprintf(fileEstimation,'%f\n',length(calls{1,nbScen}));
    for i=1:length(calls{1,nbScen})
        fprintf(fileEstimation,'%f %f %f %f %f %f %f %f %f %f %f %f\n',calls{1,nbScen}(i).time,calls{1,nbScen}(i).region,calls{1,nbScen}(i).priority,calls{1,nbScen}(i).day,calls{1,nbScen}(i).time_on_scene,calls{1,nbScen}(i).lat,calls{1,nbScen}(i).long,calls{1,nbScen}(i).time_at_hospital,calls{1,nbScen}(i).timeCleaningBase,calls{1,nbScen}(i).cleaning_needed,calls{1,nbScen}(i).hosp_needed,calls{1,nbScen}(i).ih);
    end
end
fclose(fileEstimation);



