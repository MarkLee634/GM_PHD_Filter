%GM_PHD_Initialisation
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au 

%This file initialises most of the variables that we use for the filter.

%If you want to use this GM-PHD filter for your own problem, you will need
%to replace a lot of this script with your own code.

%% Read Rosbag and fill in PX, PY
measure_bag = rosbag("jpdaf_track_corrected_2021-02-11-18-08-33.bag");
pos_bag = select(measure_bag, 'Topic', "/hummingbird1/track/bounding_box");
jpdaf_est_bag = select(measure_bag, 'Topic', "/jpdaf_tracks");
pos_array = read_array_pos(pos_bag);
jpdaf_est = read_array_est(jpdaf_est_bag);

jpdafTarget1History = [];
jpdafTarget2History = [];
jpdafTarget3History = [];
errorJPDAF_T1History = [];

%% Control Variables
%These variables control some aspect of filter function and performance
VERBOSE = 0;%Set to 1 for much more text output. This can be used to give more information about the state of the filter, but slows execution and makes the code less neat by requiring disp() statements everywhere.
KNOWN_TARGET = 1;%Set to 1 to have the first two targets already known to the filter on initialisation. Otherwise tracks for them will need to be instantiated.
PLOT_ALL_MEASUREMENTS = 0;%Set to 1 to maintain the plot of the full measurement history, including clutter and error ellipses. Set to 0 to just plot the most recent measurement.
OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS = 0;%Set to 1 to implement Vo & Ma's state extraction where targets with rounded weight greater than 1 are output multiple times. This does NOT change filter processing as extracted targets are not fed back into the filter, it is merely for display. VERBOSE must be set to 1 to see the effects of this. Personally I think this is a bit redundant but I've included it to match the Vo & Ma algorithm.
CALCULATE_OSPA_METRIC = 1; %Set to 1 to calculate the OSPA performance metric for each step. This is not essential for the GM-PHD filter but provides a rough indication of how well the filter is performing at any given timestep.
USE_EKF = 0;%Set to 1 to use extended Kalman filter. Set to 0 to use linear KF.
DRAW_VELOCITIES = 0;%Set to 1 to draw velocity arrows on the output plot.
USE_REAL_DATA = 1;
DATA_SIZE = size(pos_array(1).x, 1);

%Target initialisation: when we add a new target, we can use a two-step
%initialisation where every target is added after two observations, so it
%is added with a position (the second observation position) and a velocity
%(the velocity from the first to the second observation). Or we can add
%targets as static objects after one observation, and let the velocity be
%filled in by subsequent observations. The big difference is the first
%prediction step - with no velocity initially they will be expected to stay
%in the same location, as opposed to move in a straight line.
%We can add both types of target, but using a two-step initialisation: so a
%after the second observation, a static target is added and a target with
%velocity.
addVelocityForNewTargets = 1;
addStaticNewTargets = 0;
%When recalculating the weights of targets, we consider the match between
%observed and predicted position. We could use the observed and predicted
%velocity as well, but I haven't gotten around to implementing that. So for
%reweighting, we only use the first two values of the state.
weightDataRange = 1:2;
stateDataRange = 1:4;
temp = [];
%Filtering velocity: We throw away observed targets with velocities greater
%than this.
if USE_REAL_DATA == 0
    MAX_V = 100;
else
    MAX_V = 400;
end


%% Utility functions
calculateDataRange2 = @(j) (2*j-1):(2*j);%Used to calculate the indices of two-dimensional target j in a long list of two-dimensional targets
calculateDataRange4 = @(j) (4*(j-1)+1):(4*j);%Used to calculate the indices of four-dimensional target j in a long list of four-dimensional targets

%% Data structures and variables
%Step 0: Initialisation
k = 0; %Time step

%Step 1: Prediction for birthed targets 
numBirthedTargets = 0;
numSpawnedTargets = 0;
m_birth_before_prediction = [];%We store the birth/spawn position from before the prediction. This is used in augmenting the measurement vector to calculate velocity, for the update. 

%Step 2: Prediction for existing targets
%We store the existing position from before the prediction. This is used in augmenting measurement vector to calculate velocity, for the update.
mk_k_minus_1_before_prediction = [];%Used in augmenting measurement vector to calculate velocity, for update.
numTargets_Jk_k_minus_1 = NUMBER_OF_DRONES; %0;%Number of targets given previous. J_k|k-1. Set at end of Step 2 (prediction of existing targets)
prob_survival = 1; %0.99; %Probability of target survival. Used in GM_PHD_Predict_Existing for weight calculation

%These are used for storing the state of existing targets after prediction.
wk_k_minus_1 = [];%Weights of gaussians, previous, predicted. w_k|k-1.
mk_k_minus_1 = []; %Means of gaussians, previous, predicted. m_k|k-1
Pk_k_minus_1 = []; %Covariances of gaussians, previous, predicted. P_k|k-1

%Step 3: Construct update component.
%There is nothing we need to initialise here because the arrays are initialised
%each iteration. (This is actually true for a lot of the arrays being initialised in this script arrays, but not
%all, and rather than try to separate them it's easiest to initialise nearly all
%of them here).

%Step 4: Update
%The posterior weight, mean and covariance after the update
w_k = [];%Weights of gaussians, updated. w_k|k
m_k = [];%Means of gaussians, updated. m_k|k
P_k = [];%covariances of gaussians, updated. P_k|k
score_k = []; %Score of target
numTargets_Jk = NUMBER_OF_DRONES;%Number of targets after update. J_k. Set at end of step 4 (update)

%Step 5: Prune
numTargets_J_pruned = 0;%Number of targets after pruning
numTargets_Jk_minus_1 = NUMBER_OF_DRONES; %Number of targets, previous. J_k-1. Set in end of GM_PHD_Prune
%Merge and prune constants
%These numbers have a HUGE impact on the performance and unfortunately need to be
%manually tuned if we make changes.
%These particular values come from Vo&Ma
T = 10^-5;%Weight threshold. Value the weight needs to be above to be considered a target rather than be deleted immediately.
mergeThresholdU = 4; %Merge threshold. Points with Mahalanobis distance of less than this between them will be merged.
weightThresholdToBeExtracted = 0.5;%Value the weight needs to be above to be considered a 'real' target.

%The previous iteration's mean/weight/covariance. Set in GM_PHD_Prune
%after pruning. Used as input for the prediction step. 
wk_minus_1 = []; %Weights from previous iteration
mk_minus_1 = []; %Means from previous iteration
Pk_minus_1 = []; %Covariances from previous iteration

%Step 6: Estimate/extract states
%We store the history of all points X_k (extracted states) for plotting
%purposes. This is updated in the end of GM_PHD_Estimate
X_k_history = [];

%Step 7: Create birthed/spawned targets to append to list next iteration
%These are set in GM_PHD_Create_Birth
%Births occur at fixed positions; spawns occur at existing targets
w_birth = [];%New births' weights
m_birth = [];%New births' means
P_birth = [];%New births' covariances
w_spawn = [];%New spawns' weights
m_spawn = [];%New spawns' means
P_spawn = [];%New spawns' covariances

%Step Sim: Generate simulated measurements
%Detected clutter is a Poisson RFS
xrange = [0 400];%X range of measurements
yrange = [0 300];%Y range of measurements
V = 4 * 10^6; %Volume of surveillance region
lambda_c = 12.5 * 10^-6; %average clutter returns per unit volume (50 clutter returns over the region)
clutter_intensity = @(z_cartesian) lambda_c * V * unifpdf_2d(xrange, yrange, z_cartesian);%Generate clutter function. There are caveats to its use for clutter outside of xrange or yrange - see the comments in unifpdf_2d.m

errorT1History = [];

%Step Metric: Calculate performance metric
%The order_p and cutoff_c control filter performance; read the paper by
%Schuhmacher et al to understand these in greater detail. The parameters
%given here by default are not particularly well tuned for this problem but
%they work alright.
order_p = 1;%The order determines how punishing the metric calculation is to larger errors; as p increases, outliers are more heavily penailsed
cutoff_c = 200;%Cutoff determines the maximum error for a point.
metric_history = [];%History of the OSPA performance metric

%% MODELS for prediction and observation
%Prediction models - used in steps 1 & 2 for prediction
I2 = eye(2);%2x2 identify matrix, used to construct matrices
Z2 = zeros(2);%2x2 zero matrix, used to construct matrices
dt = 1; %One-second sampling period
F = [ [I2, dt*I2]; [Z2 I2] ];%State transition matrix (motion model)
sigma_v = 5; %Standard deviation of process noise is 5 m/(s^2)
Q = sigma_v^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; %Process noise covariance, given in Vo&Ma.

%% Birth model. This is a Poisson random finite set.
%We only use the first two elements as the initial velocities are unknown
%and it is not specified by Vo&Ma whether they are used (I am fairly sure
%they aren't, or they would have said otherwise)
%Birth and spawn models
birth_mean1 = [pos_array(1).x(1), pos_array(1).y(1), 0, 0]';%Used in birth_intensity function
birth_mean2 = [pos_array(2).x(1), pos_array(2).y(1), 0, 0]';%Used in birth_intensity function
birth_mean3 = [pos_array(3).x(1), pos_array(3).y(1), 0, 0]';%Used in birth_intensity function

covariance_birth = diag([10, 10, 5, 5]');%Used in birth_intensity function
covariance_spawn = diag([10, 10, 5, 5]');%Used in spawn_intensity function
covariance_spawn = max(covariance_spawn, 10^-6);%Used in spawn_intensity function

birth_intensity = @(x) (0.1 * mvnpdf(x(1:2)', birth_mean1(1:2)', covariance_birth(1:2,1:2)) + 0.1 * mvnpdf(x(1:2)', birth_mean2(1:2)', covariance_birth(1:2,1:2)) + (0.1 * mvnpdf(x(1:2)', birth_mean3(1:2)', covariance_birth(1:2,1:2))) );%Generate birth weight. This only takes into account the position, not the velocity, as Vo&Ma don't say if they use velocity and I assume that they don't. Taken from page 8 of their paper.
spawn_intensity = @(x, targetState) 0.05 * mvnpdf(x, targetState, covariance_spawn);%Spawn weight, from page 8 of Vo&Ma. 

prob_detection = 1; % 0.98; %Probability of target detection. Used in recalculating weights in GM_PHD_Update

%Detection models for the linear Kalman filter. The extended Kalman filter
%uses different models.
if(USE_EKF == 0)
    H = [I2, Z2];%Observation matrix for position. Not used, but if you wanted to cut back to just tracking position, might be useful.
    H2 = eye(4);%Observation matrix for position and velocity. This is the one we actually use, in GM_PHD_Construct_Update_Components
    sigma_r = 5; %Standard deviation of measurement noise is 10m. Used in creating R matrix (below)
    R = sigma_r^2 * I2;%Sensor noise covariance. used in R2 (below)
    R2 = [ [R, Z2]; [Z2, 2*R] ];%Measurement covariance, expanded to both position & velocity. Used in GM_PHD_Construct_Update_Components. NOTE: This assumes that speed measurements have the same covariance as position measurements. I have no mathematical justification for this.
end


%% read bag pos data 1,2,3
function pos_array = read_array_pos(array)

DOWN_SAMPLE = 3;

robot_posStructs = readMessages(array,'DataFormat','struct');
x = cellfun(@(m) double(m.Poses(1).Position.X),robot_posStructs);
y = cellfun(@(m) double(m.Poses(1).Position.Y),robot_posStructs);
z = cellfun(@(m) double(m.Poses(1).Position.Z),robot_posStructs);
odom_time = cellfun(@(m) double(double(m.Header.Stamp.Sec)+double(m.Header.Stamp.Nsec)*10e-10),robot_posStructs);
pos_arrayOne.x = x(1:DOWN_SAMPLE:end); %x;
pos_arrayOne.y = y(1:DOWN_SAMPLE:end); %y;
pos_arrayOne.z = z(1:DOWN_SAMPLE:end); %z;
pos_arrayOne.time = odom_time;

x2 = cellfun(@(m) double(m.Poses(2).Position.X),robot_posStructs);
y2 = cellfun(@(m) double(m.Poses(2).Position.Y),robot_posStructs);
z2 = cellfun(@(m) double(m.Poses(2).Position.Z),robot_posStructs);
odom_time2 = cellfun(@(m) double(double(m.Header.Stamp.Sec)+double(m.Header.Stamp.Nsec)*10e-10),robot_posStructs);
pos_arrayTwo.x = x2(1:DOWN_SAMPLE:end); %x2;
pos_arrayTwo.y = y2(1:DOWN_SAMPLE:end); %y2;
pos_arrayTwo.z = z2(1:DOWN_SAMPLE:end); %z2;
pos_arrayTwo.time = odom_time2;

x3 = cellfun(@(m) double(m.Poses(3).Position.X),robot_posStructs);
y3 = cellfun(@(m) double(m.Poses(3).Position.Y),robot_posStructs);
z3 = cellfun(@(m) double(m.Poses(3).Position.Z),robot_posStructs);
odom_time3 = cellfun(@(m) double(double(m.Header.Stamp.Sec)+double(m.Header.Stamp.Nsec)*10e-10),robot_posStructs);
pos_array3.x = x3(1:DOWN_SAMPLE:end); %x3;
pos_array3.y = y3(1:DOWN_SAMPLE:end); %y3;
pos_array3.z = z3(1:DOWN_SAMPLE:end); %z3;
pos_array3.time = odom_time3;

pos_array = [pos_arrayOne pos_arrayTwo pos_array3];
end

%% read jpdaf est
function est_array = read_array_est(array)

DOWN_SAMPLE = 3;

jpdaf_in = readMessages(array,'DataFormat','struct');
p1x = cellfun(@(m) double(m.Tracks(1).X),jpdaf_in);
p1y = cellfun(@(m) double(m.Tracks(1).Y),jpdaf_in);
p1.x = p1x(1:DOWN_SAMPLE:end);
p1.y = p1y(1:DOWN_SAMPLE:end);

p2x = cellfun(@(m) double(m.Tracks(2).X),jpdaf_in);
p2y = cellfun(@(m) double(m.Tracks(2).Y),jpdaf_in);
p2.x = p2x(1:DOWN_SAMPLE:end);
p2.y = p2y(1:DOWN_SAMPLE:end);

p3x = cellfun(@(m) double(m.Tracks(3).X),jpdaf_in);
p3y = cellfun(@(m) double(m.Tracks(3).Y),jpdaf_in);
p3.x = p3x(1:DOWN_SAMPLE:end);
p3.y = p3y(1:DOWN_SAMPLE:end);

est_array = [p1 p2 p3];

end