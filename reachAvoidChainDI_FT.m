total_time=tic;
%%%%%%%%%%%%%%%%%%%%%%%%%% Program begins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
concatenated_disturbance_mean=kron(ones(last_time_step,1),mean_vector);
concatenated_disturbance_sigma=kron(eye(last_time_step),sigma_matrix);

%% Reach Avoid Tube
% skip t=0; create bounds for K_1,K_2,...K_{T-1} and T
if state_dimension<4
    reachAvoidTubeLB=[kron(ones(last_time_step-1,1),safe_set_lower_bounds);target_set_LB_DP'];%-repmat(xinc,state_dimension*last_time_step,1);
    reachAvoidTubeUB=[kron(ones(last_time_step-1,1),safe_set_upper_bounds);target_set_UB_DP'];%+repmat(xinc,state_dimension*last_time_step,1);
else
    reachAvoidTubeLB=[kron(ones(last_time_step-1,1),safe_set_lower_bounds);target_set_lower_bounds];%-repmat(xinc,state_dimension*last_time_step,1);
    reachAvoidTubeUB=[kron(ones(last_time_step-1,1),safe_set_upper_bounds);target_set_upper_bounds];%+repmat(xinc,state_dimension*last_time_step,1);
end

%% concatenated_A matrix creation
concatenated_A_matrix=eye(state_dimension);
for rowNumber=1:last_time_step
    concatenated_A_matrix=[concatenated_A_matrix;system_matrix*concatenated_A_matrix(end-state_dimension+1:end,:)];
end

%% H creation
H_matrix=zeros(state_dimension*(last_time_step+1),input_dimension*last_time_step);
for rowNumber=1:last_time_step
    % Construct H_matrix block row wise --- n*(mT)
    H_matrix_temp=zeros(state_dimension,input_dimension*last_time_step);
    % What is the power of A in the first column? Ignoring the first block
    % row of zeros, it is the row number-1. (Ignoring is done later)
    maximum_exponent_for_system_matrix=rowNumber-1;
    % Construct the H_matrix block row of interest
    for system_matrix_exponent=maximum_exponent_for_system_matrix:-1:0
        column_left_indx=input_dimension*(maximum_exponent_for_system_matrix-system_matrix_exponent);
        H_matrix_temp(:,column_left_indx+1:column_left_indx+input_dimension)=system_matrix^system_matrix_exponent*input_matrix;
    end
    % The indices in the LHS ensures first row block is skipped
    H_matrix( rowNumber*state_dimension + 1 : rowNumber*state_dimension + state_dimension,:)=H_matrix_temp;    
end

%% G creation
G_matrix=zeros(state_dimension*(last_time_step+1),disturbance_dimension*last_time_step);
for rowNumber=1:last_time_step
    % Construct H_matrix block row wise --- n*(mT)
    G_matrix_temp=zeros(state_dimension,disturbance_dimension*last_time_step);
    % What is the power of A in the first column? Ignoring the first block
    % row of zeros, it is the row number-1. (Ignoring is done later)
    maximum_exponent_for_system_matrix=rowNumber-1;
    % Construct the H_matrix block row of interest
    for system_matrix_exponent=maximum_exponent_for_system_matrix:-1:0
        column_left_indx=disturbance_dimension*(maximum_exponent_for_system_matrix-system_matrix_exponent);
        G_matrix_temp(:,column_left_indx+1:column_left_indx+disturbance_dimension)=system_matrix^system_matrix_exponent*disturbance_matrix;
    end
    % The indices in the LHS ensures first row block is skipped
    G_matrix( rowNumber*state_dimension + 1 : rowNumber*state_dimension + state_dimension,:)=G_matrix_temp;    
end

%% Stochastics of G*W
concatenated_state_mean_without_input_without_state=G_matrix*concatenated_disturbance_mean;
concatenated_state_sigma=G_matrix*concatenated_disturbance_sigma*G_matrix';

%% Open-loop policy space
A_inequalities=[-eye(input_dimension*last_time_step);
                 eye(input_dimension*last_time_step)];
b_inequalities=[-umin*ones(last_time_step,1);
                 umax*ones(last_time_step,1)];

%options=optimoptions('fmincon','Display','iter','PlotFcns',@optimplotfval,'TolFun',1e-2);        
%options=optimoptions('fmincon','Display','off');%,'TolFun',1e-2); 
options=optimoptions('fmincon','Display',displayString_fmincon);   
maxDiff=0;
safe_set_center=(safe_set_lower_bounds+safe_set_upper_bounds)/2;
target_set_center=(target_set_lower_bounds+target_set_upper_bounds)/2;
chebyshevCenterRATube=[kron(ones(last_time_step-1,1),safe_set_center);target_set_center];
H_matrix_without_initial_state=H_matrix(state_dimension+1:end,:);
point_timer_FM=zeros(no_of_testing_initial_states,1);
point_timer_PS=zeros(no_of_testing_initial_states,1);
testing_RA_OP=zeros(no_of_testing_initial_states,1);
testing_RA_OP_PS=zeros(no_of_testing_initial_states,1);
testing_consv=zeros(no_of_testing_initial_states,1);
testing_consv_PS=zeros(no_of_testing_initial_states,1);
for testing_initial_state_indx=1:no_of_testing_initial_states
    initial_state=testing_initial_states(:,testing_initial_state_indx);
    X0=concatenated_A_matrix*initial_state;
    concatenated_state_mean_without_input=concatenated_state_mean_without_input_without_state+X0;
    % extract x_1,x_2,...x_N
    concatenated_state_mean_without_input_initial_state=concatenated_state_mean_without_input(state_dimension+1:end);
            
    %% U estimated to keep the mean vector close to origin (Chebyshev center of SafeTargetTube)
%     fun = @(U) norm(concatenated_state_mean_without_input_initial_state+H_matrix_without_initial_state*U-chebyshevCenterRATube,1);     
%     initial_guess_for_fmincon=fmincon(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],options);    
    fun = @(U)-log(RAprob(U,concatenated_state_mean_without_input_initial_state,concatenated_state_sigma,H_matrix_without_initial_state,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,myeps));
    if usePatternSearch==1
        point_timer_tic_PS=tic;
        PSoptions = psoptimset('Display',displayString_fmincon,'CompletePoll','on','PollMethod','GSSPositiveBasisNp1');
        [U_star_PS,RAprob_star_OLF_log_PS] = patternsearch(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],PSoptions);
        testing_RA_OP_PS(testing_initial_state_indx)=exp(-RAprob_star_OLF_log_PS);
        point_timer_PS(testing_initial_state_indx)=toc(point_timer_tic_PS);
    end
    if useFmincon==1
        point_timer_tic_FM=tic;
        [U_star,RAprob_star_OLF_log] = fmincon(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],options);                
        testing_RA_OP(testing_initial_state_indx)=exp(-RAprob_star_OLF_log);
        point_timer_FM(testing_initial_state_indx)=toc(point_timer_tic_FM);
    end
    if state_dimension<=3
        testing_consv(testing_initial_state_indx)=testing_RA_DP(testing_initial_state_indx)-testing_RA_OP(testing_initial_state_indx);
        testing_consv_PS(testing_initial_state_indx)=testing_RA_DP(testing_initial_state_indx)-testing_RA_OP_PS(testing_initial_state_indx);
%         fprintf('%d | Approx: %1.2e | True: %1.2e | Diff: %1.2e\n',testing_initial_state_indx,testing_RA_OP(testing_initial_state_indx),testing_RA_DP(testing_initial_state_indx),testing_consv(testing_initial_state_indx));
        fprintf('%d | True: %1.2e | Approx: %1.2e | Diff: %1.2e | Approx PS: %1.2e | Diff PS:%1.2e\n',testing_initial_state_indx,testing_RA_DP(testing_initial_state_indx),testing_RA_OP(testing_initial_state_indx),testing_consv(testing_initial_state_indx),testing_RA_OP_PS(testing_initial_state_indx),testing_consv_PS(testing_initial_state_indx));    
    else
        fprintf('%d | Approx: %1.2e in %1.4e seconds | Approx PS: %1.2e in %1.4e seconds\n',testing_initial_state_indx,testing_RA_OP(testing_initial_state_indx),point_timer_FM(testing_initial_state_indx),testing_RA_OP_PS(testing_initial_state_indx),point_timer_PS(testing_initial_state_indx));            
    end
end
testing_elapsed_time_FT=toc(total_time);
fprintf('Elapsed Time for %d points: %1.2e\n',no_of_testing_initial_states,testing_elapsed_time_FT);
eval(sprintf('elapsed_time_FT{%d}=testing_elapsed_time_FT;',state_dimension));
eval(sprintf('elapsed_time_FT_FM{%d}=sum(point_timer_FM);',state_dimension));
eval(sprintf('elapsed_time_FT_PS{%d}=sum(point_timer_PS);',state_dimension));
