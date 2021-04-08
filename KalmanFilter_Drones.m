
clear all;
close all;
clc;

NUMBER_OF_DRONES = 1;

%Step 0: Initialisation
GM_PHD_Initialisation_Drones;
%Prediction models - used in steps 1 & 2 for prediction
% I2 = eye(2);%2x2 identify matrix, used to construct matrices
% Z2 = zeros(2);%2x2 zero matrix, used to construct matrices
% dt = 1; %One-second sampling period
% F = [ [I2, dt*I2]; [Z2 I2] ];%State transition matrix (motion model)
% sigma_v = 5; %Standard deviation of process noise is 5 m/(s^2)
% Q = sigma_v^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; %Process noise covariance, given in Vo&Ma.
H = eye(4); %transformation matrix for cov
R = eye(4)*5; %error (noise) in observation

GM_PHD_Simulate_Initialise_Drones;

init_state = [pos_array(1).x(1), pos_array(1).y(1) 0 0]';
init_cov = [10 0 0 0; 0 10 0 0; 0 0 5 0; 0 0 0 5];

xk_k = init_state;
Pk_k = init_cov;


K_history = [];
Xk_history = [];
Xk_history_full = [];
xkp_k_history = [];
error_history = [];
zk_history = [];
Pk_k_history = [];

%Main loop
while (k < endTime)%k = timestep
   k = k + 1;
   GM_PHD_Simulate_Measurements_Drones;
   z_noise = [((rand(1)-1/2))*10 (rand(1)-1/2)*10]'; %random to give noise to perfect sim mesurement data

   
   %% predict 
   % predict state and cov
   xkp_k = F*xk_k; %+G*u (dont need since no control input)
   Pkp_k = Q + F*Pk_k*F';
   
   %% update
   %kalman gain
   Sk = H*Pkp_k*H' + R;
   K = Pkp_k*H'/Sk;
   
   %measurement
   zk = Z(:,1) + z_noise;
   zk = [zk; 0 ;0];
   
   %correct state prediction
   Xk  = xkp_k + K*(zk - xkp_k);
   %correct cov prediction
   Pk = Pkp_k - K*Sk*K';
   
   %print error
%    error = Z(:,1) - Xk(1:2)

    Pk_k_history = [Pk_k_history Pk_k];

    zk_history = [zk_history zk];
    K_history = [K_history K];
    xkp_k_history = [xkp_k_history xkp_k];
    Xk_history = [Xk_history Xk(1:2)];
    Xk_history_full = [Xk_history_full Xk];
    error = Z(:,1) - Xk(1:2);
    error = sqrt(power(Z(1,1) - Xk(1),2) + power(Z(2,1) - Xk(2),2))  ;
    error_history = [error_history error];
    figure(3)
    plot(error_history)
    
    
    %% plot
%     figure(1);
%     clf; % clear uncertainty. estimation stored in history
%     hold on;
%     axis([100 300 80 200]);
%     xlim([100 300]);
%     ylim([80 200]);
%  
%     %plot(Z(1,:), Z(2,:), 'xk');
%     %Plot target 1 track position as red dots
%     plot(Xk_history(1,1:end), Xk_history(2,1:end), '.k');
%     %Plot target 2 track position as blue dots
%     plot(zk_history(1,1:end), zk_history(2,1:end), '.r');
%     legend('estimate','measured w noise')
%     
%     
%     
%     
%     %% plot
%     %Individual X and Y components of measurements
%     figure(2);
%     hold on;
%     subplot(2,1,1);
%     plot(k, Xk_history(1,k), '.k');%X corrected estimate with kalman filter
%     plot(k, zk_history(1,k), 'or');%X measured (with random noise)
%     plot(k, xkp_k_history(1,k), 'xb');%X predicted
%     legend('KF estimate','measured', 'predict')
%     
%     subplot(2,1,2);
%     plot(k, Xk_history(2,k), '.k');%X coord of clutter measurements
%     plot(k, zk_history(2,k), 'or');%X coord of true measurement
%     plot(k, xkp_k_history(2,k), 'xb');%X coord of true measurement
%     legend('KF estimate','measured', 'predict')

   
   %update k-1 terms
   Pk_k = Pk;
   xk_k = Xk;
   
   
   


end

