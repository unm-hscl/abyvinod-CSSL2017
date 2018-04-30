tic

usePatternSearch=1;
useFmincon=1;
%%%%%%%%%%%%%%%%%%%%%%%%%% Program begins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
concatenated_disturbance_mean=kron(ones(last_time_step,1),mean_vector);
concatenated_disturbance_sigma=kron(eye(last_time_step),sigma_matrix);

%% Reach Avoid Tube
% skip t=0; create bounds for K_1,K_2,...K_{T-1} and T
reachAvoidTubeLB=[kron(ones(last_time_step-1,1),safe_set_lower_bounds);target_set_LB_DP];%-repmat(xinc,state_dimension*last_time_step,1);
reachAvoidTubeUB=[kron(ones(last_time_step-1,1),safe_set_upper_bounds);target_set_UB_DP];%+repmat(xinc,state_dimension*last_time_step,1);

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

RA_OP=zeros(length(x1), length(x2));
RA_OP_PS=zeros(length(x1), length(x2));
consv=zeros(length(x1), length(x2));
consv_PS=zeros(length(x1), length(x2));
%options=optimoptions('fmincon','Display','iter','PlotFcns',@optimplotfval,'TolFun',1e-2);        
%options=optimoptions('fmincon','Display','off');%,'TolFun',1e-2); 
options=optimoptions('fmincon','Display','off');   
maxDiff=0;
minDiff=1;
maxDiff_PS=0;
minDiff_PS=1;
safe_set_center=(safe_set_lower_bounds+safe_set_upper_bounds)/2;
target_set_center=(target_set_lower_bounds+target_set_upper_bounds)/2;
chebyshevCenterRATube=[kron(ones(last_time_step-1,1),safe_set_center);target_set_center];
H_matrix_without_initial_state=H_matrix(state_dimension+1:end,:);
for i = 1:length(x1)
        for j = 1:length(x2)
            initial_state=[x1(i);x2(j)];
            X0=concatenated_A_matrix*initial_state;
            concatenated_state_mean_without_input=concatenated_state_mean_without_input_without_state+X0;
            % extract x_1,x_2,...x_N
            concatenated_state_mean_without_input_initial_state=concatenated_state_mean_without_input(state_dimension+1:end);
            
            
            %% U estimated to keep the mean vector close to origin (Chebyshev center of SafeTargetTube)
%             fun = @(U) norm(concatenated_state_mean_without_input_initial_state+H_matrix_without_initial_state*U-chebyshevCenterRATube,1);     
%             initial_guess_for_fmincon=fmincon(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],options);    
%             fun = @(U)-log(RAprob(U,concatenated_state_mean_without_input_initial_state,concatenated_state_sigma,H_matrix_without_initial_state,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,true));
            fun = @(U)-log(RAprob(U,concatenated_state_mean_without_input_initial_state,concatenated_state_sigma,H_matrix_without_initial_state,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,myeps));
            if usePatternSearch==1
                PSoptions = psoptimset('Display','off');
                [U_star_PS,RAprob_star_OLF_log_PS] = patternsearch(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],PSoptions);
                RA_OP_PS(i,j)=exp(-RAprob_star_OLF_log_PS);
                consv_PS(i,j)=Vtsquare(i,j)-RA_OP_PS(i,j);
                if maxDiff_PS<consv_PS(i,j)
                    maxDiff_PS=consv_PS(i,j);
                end
                if minDiff_PS>consv_PS(i,j)
                    minDiff_PS=consv_PS(i,j);
                end
            end
            if useFmincon==1
                [U_star,RAprob_star_OLF_log] = fmincon(fun,zeros(input_dimension*last_time_step,1),A_inequalities,b_inequalities,[],[],[],[],[],options);                
                RA_OP(i,j)=exp(-RAprob_star_OLF_log);
                consv(i,j)=Vtsquare(i,j)-RA_OP(i,j);
                if maxDiff<consv(i,j)
                    maxDiff=consv(i,j);
                end
                if minDiff>consv(i,j)
                    minDiff=consv(i,j);
                end
            else
                RA_OP(i,j)=RA_OP_PS(i,j);
                consv(i,j)=consv_PS(i,j);
                maxDiff=maxDiff_PS;
                minDiff=minDiff_PS;
            end
            if useFmincon==1 && usePatternSearch==1
                fprintf('(%1.2f,%1.2f) | Approx: %1.2e | True: %1.2e | Diff: %1.2e | MaxDiff: %1.2e | MinDiff: %1.2e | PS -> Diff: %1.2e | MaxDiff: %1.2e | MinDiff: %1.2e\n',x1(i),x2(j),RA_OP(i,j),Vtsquare(i,j),consv(i,j),maxDiff,minDiff,consv_PS(i,j),maxDiff_PS,minDiff_PS);
            else
                fprintf('(%1.2f,%1.2f) | Approx: %1.2e | True: %1.2e | Diff: %1.2e | MaxDiff: %1.2e | MinDiff: %1.2e\n',x1(i),x2(j),RA_OP(i,j),Vtsquare(i,j),consv(i,j),maxDiff,minDiff);
            end
        end
end
elapsed_time_FT=toc;
fprintf('Elapsed Time: %1.2f\n',elapsed_time_FT)
%             [U_star,RA_star_OLF_log]=particleswarm(fun,last_time_step,umin*ones(last_time_step,1),umax*ones(last_time_step,1));
%             fun = @(U)-RAprob(U,concatenated_state_mean_without_input,concatenated_state_sigma,H_matrix,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,true);
%             RA_OP(i,j)=-RAprob_star_OLF_log;
            %% Fixed zero input policy
%             RA_OP(i,j)=RAprob(zeros(input_dimension*last_time_step,1),concatenated_state_mean_without_input,concatenated_state_sigma,H_matrix,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,true);
