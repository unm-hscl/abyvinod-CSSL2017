clear
clc
close all

%% System dynamics
time_horizon=10;
last_time_step=time_horizon;    
sampling_time=0.1;
state_dimension=2;
system_matrix=[1,sampling_time;
               0,1];
% Input
input_dimension=1;
input_matrix=[sampling_time^2/2;
              sampling_time];
umin=-1;
umax=1;
uinc=0.1;

% Disturbance
disturbance_dimension=2;
distmin=-0.5;
distmax=0.5;
distinc=0.01;
disturbance_matrix=eye(2);
mean_vector=[0;
             0];
sigma_matrix=0.01*eye(2);        


% Constraint sets
safe_corner=1;
target_corner=safe_corner/2;
safe_set_lower_bounds=safe_corner*[-1;
                                   -1];
safe_set_upper_bounds=safe_corner*[1;
                                   1];
% target_set_lower_bounds=target_corner*[-1;
%                                        0];
% target_set_upper_bounds=target_corner*[0;
%                                        1];
target_set_lower_bounds=target_corner*[-1;
                                       -1];
target_set_upper_bounds=target_corner*[1;
                                       1];
                                   
% State space
x1min = -safe_corner;     %set min position
x2min = x1min;  %set min velocity
x1max = safe_corner;      %set max position
x2max = x1max;  %set max velocity
xinc = 0.05;    %set discretization step size
x1 = x1min:xinc:x1max;      %Possible initial distances
x2 = x2min:xinc:x2max;      %Possible initial velocities

%% Target set
target_set_LB_indx=round((target_set_lower_bounds-[x1min;x2min])./xinc)+ones(2,1);
target_set_UB_indx=round((target_set_upper_bounds-[x1min;x2min])./xinc)+ones(2,1);
target_set_LB_DP=[x1(target_set_LB_indx(1));x2(target_set_LB_indx(2))];
target_set_UB_DP=[x1(target_set_UB_indx(1));x2(target_set_UB_indx(2))];

myeps=1e-2;

%% Reach Avoid using DP
% disp('Loading existing Vtsquare computation to save 60 seconds in 0.1 grid spacing and 227 seconds in 0.05 spacing');
% load('D:\Abraham\LoboDrive\MatFiles\ConsvViab\Figure2_0x05\Vtsquare_comparison')
% load('D:\LoboDrive\MatFiles\ConsvViab\Figure2_0x05\Vtsquare_comparison')
reachAvoidDI_DP

figure(1)
clf
colorbar;
hold on;
surf(wdisc1D,wdisc1D,reshape(pdisc,length(wdisc1D),[]));
set(gca,'Fontsize',15);
title('Disturbance density approximated');
xlabel('$w_1$','interpreter','latex','Fontsize',20);
ylabel('$w_2$','interpreter','latex','Fontsize',20);
box on;
grid on;
view([0,90]);
drawnow

%% Reach Avoid using FourierTransforms
reachAvoidDI_FT

disp('Comment the line that loads the matfile in Figure 2.m and run it to get Figure 2 based on this data');