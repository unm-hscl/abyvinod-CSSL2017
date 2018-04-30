clear
clc
close all

%% System dynamics
no_of_testing_initial_states=20;
time_horizon=10;
last_time_step=time_horizon;    
sampling_time=0.1;
state_dimension=2;
test_point_x=0.1;
test_point_y=0.9;

% Input
input_dimension=1;
umin=-1;
umax=1;
uinc=0.1;
% Input gridding
udisc=umin:uinc:umax;
no_of_input_points=length(udisc);

% Disturbance
distmin=-0.5;
distmax=0.5;
distinc=0.05;
sigma_square=0.01;

% Constraint sets
adjustment_for_random_initial_test_condition=5;
safe_corner=1;
target_corner=safe_corner/2;
myeps=1e-3;

% State space
xmin = -safe_corner;     %set min position
xmax = safe_corner;      %set max position
xinc = 0.5;    %set discretization step size
elapsed_time_DP={};
elapsed_time_FT={};
elapsed_time_FT_FM={};
elapsed_time_FT_PS={};
initial_states_for_testing={};
run_DP_script=1;

if run_DP_script==0
    load('D:\LoboDrive\MatFiles\ConsvViab\curseOfDim_full_10_5_100points');
end

grid_RA_DP_val={};
indx_grid_RA_DP_val=1;
compute_time_grid=[];
for xinc=[0.1,0.05,0.01,0.005]
    fprintf('Solving for state_grid_spacing=%f\n',xinc);
    x = xmin:xinc:xmax;      %Possible initial distances
    indices=[sub2ind([length(x),length(x)],find(abs(x-test_point_x)<xinc-1e-16),find(abs(x-test_point_y)<xinc-1e-16)),sub2ind([length(x),length(x)],find(abs(x+test_point_x)<xinc-1e-16),find(abs(x+test_point_y)<xinc-1e-16))];
    safe_set_lower_bounds=-safe_corner*ones(state_dimension,1);
    safe_set_upper_bounds=safe_corner*ones(state_dimension,1);
    target_set_lower_bounds=-target_corner*ones(state_dimension,1);
    target_set_upper_bounds=target_corner*ones(state_dimension,1);

    % System and input matrices 
    system_matrix=eye(state_dimension);
    input_matrix=zeros(state_dimension,1);
    disturbance_dimension=state_dimension;
    disturbance_matrix=eye(state_dimension);
    mean_vector=zeros(state_dimension,1);
    sigma_matrix=sigma_square*eye(state_dimension);        
    for i=1:state_dimension
        for j=i+1:state_dimension
            system_matrix(i,j)=sampling_time^(j-i)/factorial(j-i);
        end
        input_matrix(i,1)=sampling_time^(state_dimension-i+1)/factorial(state_dimension-i+1);
    end
    
    %% Reach Avoid using DP
    target_set_LB_indx=round((target_set_lower_bounds-repmat(xmin,state_dimension,1))./xinc)+ones(state_dimension,1);
    target_set_UB_indx=round((target_set_upper_bounds-repmat(xmin,state_dimension,1))./xinc)+ones(state_dimension,1);
    target_set_LB_DP=x(target_set_LB_indx);
    target_set_UB_DP=x(target_set_UB_indx);

    reachAvoidChainDI_DP          
    testing_RA_DP=terminal_RA_prob{state_dimension}(indices)
    grid_RA_DP_val{indx_grid_RA_DP_val}=testing_RA_DP;
    compute_time_grid(indx_grid_RA_DP_val)=elapsed_time_DP{state_dimension};
    indx_grid_RA_DP_val=indx_grid_RA_DP_val+1;    
end
%% Table I contents
no_of_testing_initial_states=length(indices);
[indx_x,indx_y]=ind2sub([length(x),length(x)],indices);
testing_initial_states=[x(indx_x)',x(indx_y)']';
reachAvoidChainDI_FT
% save('D:/LoboDrive/MatFiles/ConsvViab/curse_useful.mat','terminal_RA_prob','initial_states_for_testing','initial_states_indx_for_testing','consv','RA_DP','RA_OP','elapsed_time_DP','elapsed_time_FT')
