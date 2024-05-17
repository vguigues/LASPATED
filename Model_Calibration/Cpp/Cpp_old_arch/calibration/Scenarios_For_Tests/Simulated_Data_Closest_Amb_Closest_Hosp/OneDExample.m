%All calls of the lowest priority
calls(1).time=0;
calls(1).lat=0;
calls(1).long=0;
calls(1).time_on_scene=2;
calls(1).time_at_hospital=2;

calls(2).time=1;
calls(2).lat=-5;
calls(2).long=0;
calls(2).time_on_scene=2;
calls(2).time_at_hospital=2;

calls(3).time=4;
calls(3).lat=2.5;
calls(3).long=0;
calls(3).time_on_scene=1;
calls(3).time_at_hospital=1;

calls(4).time=5;
calls(4).lat=0;
calls(4).long=0;
calls(4).time_on_scene=1;
calls(4).time_at_hospital=1;

nb_ambulances=2;

ambulances(1).base.lat=-2;
ambulances(1).base.long=0;
ambulances(1).speed=1;
ambulances(1).hospital_destination.lat=[]; 
ambulances(1).hospital_destination.long=[]; 
ambulances(1).arrival_time_at_h_last_trip=-1;
ambulances(1).arrival_time_at_b_last_trip=0;

ambulances(2).base.lat=3;
ambulances(2).base.long=0;
ambulances(2).speed=1;
ambulances(2).hospital_destination.lat=[]; 
ambulances(2).hospital_destination.long=[]; 
ambulances(2).arrival_time_at_h_last_trip=-1;
ambulances(2).arrival_time_at_b_last_trip=0;

hospitals(1).lat=-1;
hospitals(1).long=0;
hospitals(2).lat=-6;
hospitals(2).long=0;
hospitals(3).lat=2;
hospitals(3).long=0;
