%GM_PHD_Filter
%Version 1.10, last modified 7th January 2014

clear all;
close all;
clc;

NUM_DRONES = 3;

%Step 0: Initialisation
%the linear KF is direct observations of target state, the EKF is
%range-bearing measurements. 
GM_PHD_Initialisation_drones;
GM_PHD_Simulate_Initialise_drones;


%In Vo&Ma, the targets are known at filter initialisation.
if KNOWN_TARGET == 1
    t1start = [simTarget1Start(1:2); simTarget1Vel];
    t2start = [simTarget2Start(1:2); simTarget2Vel];
    t3start = [simTarget3Start(1:2); simTarget3Vel];
    
    m_birth = [t1start, t2start, t3start];
    w_birth = [birth_intensity(t1start), birth_intensity(t2start), birth_intensity(t3start)];
    P_birth = [covariance_birth, covariance_birth, covariance_birth];
    numBirthedTargets = NUM_DRONES;
end


%% Control parameters
%prune step
mergeThresholdU = 5;
%sim data noise
noiseScaler = 0.0;       %Adjust the strength of the noise on the measurements by adjusting this. Useful for debugging.
%false positive
nClutter = 0; %Assume constant 50 clutter measurements. Since clutter is Poisson distrbuted it might be more accurate to use nClutter = poissrnd(50) if you have the required Matlab toolbox. Constant 50 clutter works well enough for simulation purposes.
%false negative

%Main loop
while (k < endTime)%k = timestep
    k = k + 1;
    s = sprintf('======ITERATION %d======', k);
    disp(s);
    
    if(k >= 25)
        stophere = 0;
    end
        
    %Step 0: Sim Generate sensor Measurements
    
    GM_PHD_Simulate_Measurements_drones;  %Linear KF measurements are simulated direct observations [X; Y] of the target positions
    
    %Step 1: Prediction for birthed/spawned targets
    GM_PHD_Predict_Birth;
    %Step 2: Prediction for existing targets
    GM_PHD_Predict_Existing;
    %Step 3: Construction of PHD update components
    GM_PHD_Construct_Update_Components;
    %Step 4: Update targets with measurements
    GM_PHD_Update;
    
    
    %Step 5: Prune targets
    GM_PHD_Prune;
    %Step 6: Estimate position of targets
    GM_PHD_Estimate
    
    %Step 7: Create birthed-targets-list to add next iteration in Step 1.
    %Not a formal part of Vo&Ma but an essential step!
    %The EKF version uses an inverse sensor model.
    
    GM_PHD_Create_Birth;
    
    
    %Step Plot: Generate graphs
    
    GM_PHD_Simulate_Plot;
    



end

