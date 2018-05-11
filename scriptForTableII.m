clear
clc
close all

fprintf(strcat('We run DP for the viability problem of a double',...
        ' integrator at various grid step-sizes [0.1,0.05,0.01,0.005].',...
        '\nDue to the structure in the problem, the true optimal value',...
        ' function is even and we check for the same.\n\n'));
disp('Generate Table 2 (Choose the following options)');
code_flow_flag = input('1-Compute DP for various grids | 2-Load the matfile\n');
if code_flow_flag == 2
    load('table2_contents.mat');
elseif code_flow_flag == 1
    grid_step_sizes=[0.1,0.05,0.01,0.005];
    
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
    initial_states_for_testing={};

    grid_RA_DP_value_cell={};
    indx_grid_RA_DP_val=1;
    compute_time_grid=[];
    for xinc=grid_step_sizes
        fprintf('Solving for state_grid_spacing=%1.3f\n',xinc);
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
        testing_RA_DP=terminal_RA_prob{state_dimension}(indices);
        grid_RA_DP_value_cell{indx_grid_RA_DP_val}=testing_RA_DP;
        compute_time_grid(indx_grid_RA_DP_val)=elapsed_time_DP{state_dimension};
        indx_grid_RA_DP_val=indx_grid_RA_DP_val+1;    
    end
else
    error('Invalid option');
end
disp('Grid sizes');
disp(grid_step_sizes)
disp('Value function tested at x and -x');
disp([grid_RA_DP_value_cell{:}])
disp('Computation time');
disp(compute_time_grid)