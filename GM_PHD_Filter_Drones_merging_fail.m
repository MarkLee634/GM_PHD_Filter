%GM_PHD_Filter
%Version 1.10, last modified 7th January 2014
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au 
%With:
%- some Kalman filter update code by Tim Bailey, taken from his website http://www-personal.acfr.usyd.edu.au/tbailey/software/
%- error_ellipse by AJ Johnson, taken from Matlab Central http://www.mathworks.com.au/matlabcentral/fileexchange/4705-errorellipse

%The GM-PHD algorithm is from Ba-Ngu Vo & Wing-Kin Ma in:
%B.-N. Vo, W.-K. Ma, "The Gaussian Mixture Probability Hypothesis Density
%Filter", IEEE Transactions on Signal Processing, Vol 54, No. 11, November 2006, pp4091-4104

%I have implemented both the linear Kalman filter and extended Kalman
%filter versions of this algorithm; switch between these by
%setting/clearing the variable USE_EKF in GM_PHD_Initialisation.

%Ba-Ngu Vo has kindly allowed me to include his implementation of the Optimal Subpattern Assignment
%(OSPA) metric proposed by D. Schuhmacher, Ba-Tuong Vo & Ba-Ngu Vo in
% Schuhmacher, D.; Ba-Tuong Vo; Ba-Ngu Vo, "A Consistent Metric for Performance Evaluation of Multi-Object Filters," Signal Processing, IEEE Transactions on , vol.56, no.8, pp.3447,3457, Aug. 2008
%This uses
%- ospa_dist by Ba-Ngu Vo, taken from http://ba-ngu.vo-au.com/vo/OSPA_for_Tracks.zip
%- Hungarian by Alex Melin (to whom I am also much obliged), also taken from http://ba-ngu.vo-au.com/vo/OSPA_for_Tracks.zip

%The OSPA metric is not essential for the functioning of the filter but
%provides a nice way of analysing performance.

%See the README.txt, the comments and the Vo & Ma paper for more
%information about what this code actually does.

clear all;
close all;
clc;

USE_DRONES = 1;
NUMBER_OF_DRONES = 3;

%Step 0: Initialisation
%The EKF version must initialise the Jacobian functions used to linearise
%the observation/prediction models. The simulated measurements are also very
%different; the linear KF is direct observations of target state, the EKF is
%range-bearing measurements. The observation and prediction covariances are
%also different.

if USE_DRONES == 0
    GM_PHD_Initialisation;
else
    GM_PHD_Initialisation_Drones;
end


if USE_EKF == 0
    if USE_DRONES == 0
        GM_PHD_Simulate_Initialise;
    else
        GM_PHD_Simulate_Initialise_Drones;
    end
end

%In Vo&Ma, the targets are known at filter initialisation.
%If we want to know about them, set KNOWN_TARGET to 1 in GM_PHD_Initialisation.
%Otherwise they should be initialised after being detected a few times.
%HOWEVER this is not guaranteed - sometimes due to noise or missed
%detections, one or both of the targets will not be tracked. This is just 
%part of the filter and the birth_intensity function.
if KNOWN_TARGET == 1
    t1start = [simTarget1Start(1:2); simTarget1Vel];
    t2start = [simTarget2Start(1:2); simTarget2Vel];
    t3start = [simTarget3Start(1:2); simTarget3Vel];
    m_birth = [t1start, t2start, t3start];
    w_birth = [birth_intensity(t1start), birth_intensity(t2start), birth_intensity(t3start)];
    P_birth = [covariance_birth, covariance_birth, covariance_birth];
    numBirthedTargets = NUMBER_OF_DRONES;
end

mk_minus_1 = m_birth;
wk_minus_1 = w_birth; 
Pk_minus_1 = P_birth;

%Main loop
while (k < endTime)%k = timestep
    k = k + 1;
    s = sprintf('======ITERATION %d======', k);
    disp(s);
        
    %Step Sim: Generate sensor Measurements
    %If you want to use this code with your own data or for a different problem,
    %replace this function with your own.
    if USE_EKF == 0
        if USE_DRONES == 0
            GM_PHD_Simulate_Measurements;  %Linear KF measurements are simulated direct observations [X; Y] of the target positions
        else
            GM_PHD_Simulate_Measurements_Drones;
        end
    end
    
    %% predict birth
    % dont do anything because m,w,p initialized spawn not needed
    
    
    %% predict existing target
    mk_k_minus_1_before_prediction = mk_minus_1;
    
    for j = 1:size(mk_minus_1,2) 
        
        %update w,m,p
        wk_minus_1(j) = prob_survival * wk_minus_1(j); %maintain
        mk_minus_1(:,j) = F * mk_minus_1(:,j);
        
        P_range = calculateDataRange4(j);
        P_i = Q + F * Pk_minus_1(:,P_range) * F';
        Pk_minus_1(:,P_range) = P_i;
        
        prevState = mk_k_minus_1_before_prediction(:,j);
        newState = mk_minus_1(:,j);
        
        s = sprintf('\t\tExisting target %d. Previously at %3.1f %3.1f, now at %3.1f %3.1f.', j, prevState(1), prevState(2), newState(1), newState(2));
        disp(s);
        
    end
    
    wk_k_minus_1 = [wk_minus_1];
    mk_k_minus_1 = [mk_minus_1];
    Pk_k_minus_1 = [Pk_minus_1];
    numTargets_Jk_k_minus_1 = size(mk_minus_1,2);  
    
    
    %% construct
    eta = [];
    S = [];
    K = [];
    P_k_k = [];
    
    for j = 1:numTargets_Jk_k_minus_1 
        m_j = mk_k_minus_1(:,j);
        eta_j = H2 * m_j;
        eta = [eta, eta_j];
        
        P_range = calculateDataRange4(j); %4x4 array

        PHt = Pk_k_minus_1(:,P_range) * H2'; %Taken from Tim Bailey's EKF code. 4x4 array

        %Calculate K via Tim Bailey's method.
        S_j = R2 + H2 * PHt;
        S = [S, S_j];
        %At this point, Tim Bailey's code makes S_j symmetric. In this case, it leads to the matrix being non-positive definite a lot of the time and chol() crashes.
        %So we won't do that. 
        SChol= chol(S_j);

        SCholInv= SChol \ eye(size(SChol)); % triangular matrix, invert via left division
        W1 = PHt * SCholInv;

        K_j = W1 * SCholInv';
        K = [K, K_j];

        P_j = Pk_k_minus_1(:,P_range) - W1*W1';%4x4 array
        P_k_k = [P_k_k, P_j];
        %End Tim Bailey's code.
        
    end
    
    %% update
    w_k = zeros(1, numTargets_Jk_k_minus_1* size(Z,2) + numTargets_Jk_k_minus_1); %1x12 (3 noise + 3 target x 3 measure)
    m_k = zeros(4, numTargets_Jk_k_minus_1* size(Z,2) + numTargets_Jk_k_minus_1); %4x12
    P_k = zeros(size(Pk_k_minus_1,1), size(Pk_k_minus_1,1) * (numTargets_Jk_k_minus_1 * size(Z, 2) + numTargets_Jk_k_minus_1) );
    
    L = 0;
    
    for zi = 1:size(Z,2) %for all observed measurements
       L = L +1;
       
       for j = 1:numTargets_Jk_k_minus_1 % for all tracked targets
           
           %get velocity
           thisZ = Z(:,zi);
           prevX = mk_k_minus_1_before_prediction(1:2,j);
           thisV = (thisZ - prevX)/dt;
           thisZ = [thisZ; thisV];
           
           thisIndex = L * numTargets_Jk_k_minus_1 + j; % 4~12
           
           S_range = calculateDataRange4(j);
           %update weight
           w_new = wk_k_minus_1(j) * mvnpdf(thisZ(1:2), eta(1:2,j), S(1:2,S_range(1:2)));
           w_k(thisIndex) = w_new;
           %update mean
           delta = thisZ - eta(:,j);
           K_range = calculateDataRange4(j);
           m_new = mk_k_minus_1(:,j) + K(:,K_range) * (delta);
           m_k(:,thisIndex) = m_new;
           %update cov
           old_P_range = calculateDataRange4(j);%Returns 4 columns
           new_P_range = 4 * L * numTargets_Jk_k_minus_1 + old_P_range;
           P_new = P_k_k(:,old_P_range);
           P_k(:,new_P_range) = P_new;
           
  
       end
       
       weight_tally = 0;
       for i = 1:numTargets_Jk_k_minus_1 %sum weight to normalize
           thisIndex = L * numTargets_Jk_k_minus_1 + i;
           weight_tally = weight_tally + w_k(thisIndex);
       end
       
       for j = 1:numTargets_Jk_k_minus_1 %recalculate weight to normalize
          
           old_weight = w_k(L * numTargets_Jk_k_minus_1 + j);
           measZ = [Z(1,zi), Z(2,zi)];
           new_weight = old_weight / (clutter_intensity(measZ) + weight_tally);%Normalise    
           w_k(L * numTargets_Jk_k_minus_1 + j) = new_weight;
           
       end
       
    end
    
    numTargets_Jk = L * numTargets_Jk_k_minus_1 + numTargets_Jk_k_minus_1; %12
    
    for j = 1:numTargets_Jk
        thisPos = m_k(:,j);
        thisW = w_k(j);
%         s = sprintf('Target %d: %3.1f %3.1f %3.1f %3.1f, Weight %3.9f', j, thisPos(1), thisPos(2), thisPos(3), thisPos(4), thisW);
%         disp(s);
    end

    
    
    %% prune
    
    
    %store w,m,p for next iteration
    I = [];
    m_bar_k = zeros(4,numTargets_Jk_k_minus_1);
    w_bar_k = zeros(1,numTargets_Jk_k_minus_1);
%     P_bar_k = zeros(4,numTargets_Jk_k_minus_1*4);
%     m_bar_k = [];
%     w_bar_k = [];
    P_bar_k = [];
    
    WEIGHT_THRESHOLD_CLUSTER = 0.004;
    MERGE_THRESHOLD_DISTANCE = 20;
    % threshold decent weights
    I = find(w_k >= WEIGHT_THRESHOLD_CLUSTER)

    
    while isempty(I) == false
        
    %Find j, which is i corresponding to highest w_k for all i in I
    highWeights = w_k(:,I);
    [maxW, j] = max(highWeights);
    j = j(1); %In case of two targets with equal weight
    %j is an index of highWeights (i.e. a position in I)
    %We want the index in w_k
    j = I(j);
    
    
       % cluster if within close distance
       L = []; %index for w_k that's good weight & close
       for iterateI = 1:size(I,2) 
           thisI = I(iterateI);
           
           %get distance from maxWeight to other goodWeights
           delta_m = m_k(:,j) - m_k(:,thisI);
           mahal_dist = sqrt(delta_m' * delta_m);
%            P_range = calculateDataRange4(thisI); 
%            mahal_dist = delta_m' * (P_k(:,P_range) \ delta_m);%Invert covariance via left division
           if(mahal_dist <= MERGE_THRESHOLD_DISTANCE)
               L = [L, thisI];
               sprintf('adding %d distance %f', thisI,mahal_dist)
           end
       end %got index of close distance weights in L
       
       w_bar_k_l = sum(w_k(L));
       % states = w * m (TO DO)
       %The new mean is the weighted average of the merged means
       m_bar_k_l = 0;
       for i = 1:size(L,2)
            thisI = L(i);
            m_bar_k_l = m_bar_k_l + 1 / w_bar_k_l *  (w_k(thisI) * m_k(:,thisI));
       end
        
       %store what association this grouping was for state N
       thisIndexAssociation = floor((L(1)-1)/numTargets_Jk_minus_1);

       
        P_val = [];
        for i=1:size(I,2)
            thisI = I(i);
            P_range = calculateDataRange4(thisI);
            P_val = [P_val, P_k(:,P_range)];
        end
        
        %Now delete the elements in L from I
        for i = 1:length(L)
            iToRemove = find(I == L(i));
            I(iToRemove) = [];
        end
        
        %append to the correct target
        m_bar_k(:,thisIndexAssociation) = m_bar_k_l
        w_bar_k(:,thisIndexAssociation) = w_bar_k_l;
        
%         m_bar_k = [m_bar_k, m_bar_k_l]
%         w_bar_k = [w_bar_k, w_bar_k_l]
        
        P_bar_k = [P_bar_k, P_val];
    
   
    end
    
    numTargets_Jk_minus_1 = size(w_bar_k,2);
    
    % update states
    Pk_minus_1 = P_bar_k;
    wk_minus_1 = w_bar_k;
    mk_minus_1 = m_bar_k;
    
    
    %% state extraction
    X_k = [];
    X_k_w = [];
    X_k_P = [];
    
    X_k = m_bar_k;
    X_k_w = w_bar_k;
    X_k_P = P_bar_k;
    
    if k > 10 
       
        highWeightIndex = w_bar_k > 0.05;
        X_k = m_bar_k(:,highWeightIndex);
        X_k_w = w_bar_k(:,highWeightIndex);
        X_k_P = P_bar_k;
        
    end
    
    
    
    %Store history for plotting.
    X_k_history = [X_k_history, X_k];

    %% plot
    figure(1);
    clf; % clear uncertainty. estimation stored in history
    hold on;
    axis([0 300 0 250]);
    xlim([0 300]);
    ylim([0 250]);
 
    plot(Z(1,:), Z(2,:), 'xk');
    

    %Plot target 1 track position as red dots
    plot(X_k_history(1,1:3:end), X_k_history(2,1:3:end), '.-r');
    %Plot target 2 track position as blue dots
    plot(X_k_history(1,2:3:end), X_k_history(2,2:3:end), '.-b');
    %Plot target 3 track position as green dots
    plot(X_k_history(1,3:3:end), X_k_history(2,3:3:end), '.-g');
    
    xlabel('X position');
    ylabel('Y position');
    title('Simulated targets and measurements');
    axis square;
    
    %For extracted targets, plot latest target(s) as magenta circle
    %and draw an error ellipse to show uncertainty
    if(~isempty(X_k))
        plot(X_k(1,:), X_k(2,:), 'om');

        for c = 1: size(X_k,2)
           thisMu = X_k(1:2, c);
           covRange = calculateDataRange4(c);
           thisCov = X_k_P(:,covRange);
           thisCov = thisCov(1:2, 1:2); %We only care about position
           error_ellipse(thisCov, thisMu);
        end

    end


    %% original steps
%     %Linear KF use fixed matrices for prediction and update.
%     %Step 1: Prediction for birthed/spawned targets
%     GM_PHD_Predict_Birth;
%     %Step 2: Prediction for existing targets
%     GM_PHD_Predict_Existing;
%     %Step 3: Construction of PHD update components
%     GM_PHD_Construct_Update_Components;
%     %Step 4: Update targets with measurements
%     GM_PHD_Update;
%     %Step 5: Prune targets
%     GM_PHD_Prune;
%     %Step 6: Estimate position of targets
%     GM_PHD_Estimate
%     %Step 7: Create birthed-targets-list to add next iteration in Step 1.
%     %Not a formal part of Vo&Ma but an essential step!
%     GM_PHD_Create_Birth; 
%     %Step Metric: Calculate performance metric
%     GM_PHD_Calculate_Performance_Metric;   
%     %Step Plot: Generate graphs
%     GM_PHD_Simulate_Plot;
 %% 
 
    
    if(VERBOSE == true)
        pause;%Pause to allow reading of the text
    end

end

