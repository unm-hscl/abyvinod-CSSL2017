%% System dynamics
no_of_testing_initial_states=20;
time_horizon=10;
last_time_step=time_horizon;    
sampling_time=0.1;
state_dimension=2;

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
safe_corner=10;
target_corner=safe_corner/2;
myeps=1e-2;
if table_1_generation==1
    target_corner=safe_corner/10*8;
    myeps=1e-3;
    state_dimension_list = 40;
else
    state_dimension_list = 1:40;
end


% State space
xmin = -safe_corner;     %set min position
xmax = safe_corner;      %set max position
xinc = 0.5;    %set discretization step size
x = xmin:xinc:xmax;      %Possible initial distances
elapsed_time_DP={};
elapsed_time_FT={};
elapsed_time_FT_FM={};
elapsed_time_FT_PS={};
RA_DP={};
RA_OP={};
RA_OP_PS={};
consv={};
consv_PS={};
initial_states_for_testing={};

for state_dimension=state_dimension_list
    fprintf('Solving for state_dimension=%d\n',state_dimension);
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
    if state_dimension<4
        reachAvoidChainDI_DP
%             figure(state_dimension*10+1)
%             clf
%             colorbar;
%             hold on;
%             surf(x,x,Vtsquare');
%             set(gca,'Fontsize',15);
%             title('ReachAvoid DP approximation');
%             xlabel('$x_1$','interpreter','latex','Fontsize',20);
%             ylabel('$x_2$','interpreter','latex','Fontsize',20);
%             box on;
%             grid on;
%             view([0,90]);
%             drawnow
        %% Reach Avoid using FourierTransforms            
        indices=round((length(x)-1)*rand(state_dimension,no_of_testing_initial_states))+1;      
        testing_initial_states=x(indices);
    else
        %% Reach Avoid using FourierTransforms            
        indices=target_set_LB_indx(1)+adjustment_for_random_initial_test_condition+round((target_set_UB_indx(1)-target_set_LB_indx(1)-2*adjustment_for_random_initial_test_condition)*rand(state_dimension,no_of_testing_initial_states));
        testing_initial_states=x(indices);
    end   

    if table_1_generation
        no_of_testing_initial_states=3;
        %indices=[repmat(round((length(x)-1)/2)+1,state_dimension,1),repmat(round((length(x)-1)/2)+4,state_dimension,1),[36;5;repmat(round((length(x)-1)/2)+5,state_dimension-2,1)]];            % For n=20
        indices=[repmat(round((length(x)-1)/2)+1,state_dimension,1),repmat(round((length(x)-1)/2)+6,state_dimension,1),[4;37;repmat([4;37],state_dimension/2-1,1)]];             % For n=40
        testing_initial_states=x(indices);
    end
%     
    %% Table I contents
%     indices=[   81
%                 82
%                 83
%                 84
%                 85
%                 86
%                 87
%                 88
%                314
%                315
%                316
%                317
%                318
%                319
%                320
%                321]';
%     no_of_testing_initial_states=length(indices);
%     testing_initial_states=x(indices);

    if state_dimension>3 | run_DP_script == 0
        reshaping_scaling_matrix=1;
        for i=2:state_dimension
            reshaping_scaling_matrix=[reshaping_scaling_matrix,length(x)^(i-1)];
        end
    end
    indx_columnwise=reshaping_scaling_matrix*(indices-[zeros(1,no_of_testing_initial_states);ones(state_dimension-1,no_of_testing_initial_states)]);
    eval(sprintf('initial_states_for_testing{%d}=testing_initial_states;',state_dimension));
    eval(sprintf('initial_states_indx_for_testing{%d}=indices;',state_dimension));    
    if state_dimension<=3
        eval(sprintf('testing_RA_DP=terminal_RA_prob{%d};',state_dimension));
        testing_RA_DP=testing_RA_DP(indx_columnwise);
    end
    reachAvoidChainDI_FT
    eval(sprintf('RA_OP{%d}=testing_RA_OP;',state_dimension));
    eval(sprintf('RA_OP_PS{%d}=testing_RA_OP_PS;',state_dimension));
    if state_dimension<4
        eval(sprintf('consv{%d}=testing_consv;',state_dimension));    
        eval(sprintf('consv_PS{%d}=testing_consv_PS;',state_dimension));    
        eval(sprintf('RA_DP{%d}=testing_RA_DP;',state_dimension));
    end
end
