%GM_PHD_Simulate_Measurements
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au 

%This file generates simulated measurement data for  the simulation
%described in Vo&Ma.
%There will be gaussian noise on the measurement and Poisson-distributed clutter
%in the environment. 

%If you want to use this PHD filter implementation for another problem, you
%will need to replace this script with another one that populates Z,
%zTrue, and simMeasurementHistory (Z is used in a lot of the update code,
%zTrue and simMeasurementHistory are used in GM_PHD_Simulate_Plot)

%Note: It is possible to get no measurements if the target is not detected
%and there is no clutter
s = sprintf('Step Sim: Simulating measurements.');
disp(s);


%Simulate target movement
if k > 1
   simTarget1Vel = [(pos_array(1).x(k)-pos_array(1).x(k-1))/dt; (pos_array(1).y(k)-pos_array(1).y(k-1))/dt];
   simTarget2Vel = [(pos_array(2).x(k)-pos_array(2).x(k-1))/dt; (pos_array(2).y(k)-pos_array(2).y(k-1))/dt];
   simTarget3Vel = [(pos_array(3).x(k)-pos_array(3).x(k-1))/dt; (pos_array(3).y(k)-pos_array(3).y(k-1))/dt];   
else
    simTarget1Vel = [0;0];
    simTarget2Vel = [0;0];
    simTarget3Vel = [0;0];
end
simTarget1State = [pos_array(1).x(k); pos_array(1).y(k); simTarget1Vel]; 
simTarget2State = [pos_array(2).x(k); pos_array(2).y(k); simTarget2Vel]; 
simTarget3State = [pos_array(3).x(k); pos_array(3).y(k); simTarget3Vel]; 

jpdafTarget1State = [jpdaf_est(1).x(k); jpdaf_est(1).y(k)];
jpdafTarget2State = [jpdaf_est(2).x(k); jpdaf_est(2).y(k)];
jpdafTarget3State = [jpdaf_est(3).x(k); jpdaf_est(3).y(k)];


%Save target movement for plotting
simTarget1History = [simTarget1History, simTarget1State];
simTarget2History = [simTarget2History, simTarget2State];
simTarget3History = [simTarget3History, simTarget3State];

jpdafTarget1History = [jpdafTarget1History, jpdafTarget1State];
jpdafTarget2History = [jpdafTarget2History, jpdafTarget2State];
jpdafTarget3History = [jpdafTarget3History, jpdafTarget3State];

Z1_true = [ simTarget1State(1);  simTarget1State(2)];
Z2_true = [ simTarget2State(1);  simTarget2State(2)];
Z3_true = [ simTarget3State(1);  simTarget3State(2)];

 measX1 = simTarget1State(1);
    measY1 = simTarget1State(2);
    measX2 = simTarget2State(1);
    measY2 = simTarget2State(2);
    measX3 = simTarget3State(1);
    measY3 = simTarget3State(2);
    

%% drop measurements

%     if rem(k,10)==0
%         
%         measX1 = simTarget1State(1);
%         measY1 = simTarget1State(2);
%         measX2 = simTarget3State(1);
%         measY2 = simTarget3State(2);
%         measX3 = 0;
%         measY3 = 0;
%         
%     end


%% mix measurement order 

% if k > 120
%     measX1 = simTarget1State(1);
%     measY1 = simTarget1State(2);
%     measX2 = simTarget2State(1);
%     measY2 = simTarget2State(2);
%     measX3 = simTarget3State(1);
%     measY3 = simTarget3State(2);
% 
% elseif k > 60
%     measX2 = simTarget1State(1);
%     measY2 = simTarget1State(2);
%     measX3 = simTarget2State(1);
%     measY3 = simTarget2State(2);
%     measX1 = simTarget3State(1);
%     measY1 = simTarget3State(2);
%     
% else
%     
%     measX1 = simTarget1State(1);
%     measY1 = simTarget1State(2);
%     measX2 = simTarget2State(1);
%     measY2 = simTarget2State(2);
%     measX3 = simTarget3State(1);
%     measY3 = simTarget3State(2);
% 
% end

   

    
%Generate true measurement
Z = [ [measX1 measX2 measX3]; [measY1 measY2 measY3] ];
zTrue = Z;%Store for plotting

clutter = [];
%Append clutter
Z = [Z, clutter];

%Store history
simMeasurementHistory{k} =  Z;