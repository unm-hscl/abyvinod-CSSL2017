function prob=RAprob(U,concatenated_state_mean_without_initial_state,concatenated_state_sigma,H_matrix_without_initial_state,reachAvoidTubeLB,reachAvoidTubeUB,state_dimension,myeps)    
    mu_vector=concatenated_state_mean_without_initial_state+H_matrix_without_initial_state*U;
    sigma_matrix=concatenated_state_sigma(state_dimension+1:end,state_dimension+1:end);
    %% QSCMVNV in a loop using the error estimate
    error_quadrature=10;
    points_base=10;
    points_power=1;
    while abs(error_quadrature)>myeps
        [prob,error_quadrature]=qscmvnv( points_base^points_power, sigma_matrix, reachAvoidTubeLB-mu_vector, eye(size(mu_vector,1),size(mu_vector,1)), reachAvoidTubeUB-mu_vector);
        prob=round(prob/myeps)*myeps;
        % If zero, then saturate it for log
        if prob<myeps
            prob=myeps;
        end
        if points_power>5
            fprintf('Exceeded 5 iterations --- Required accuracy: %1.2e | Current error: %1.2e\n',myeps,error_quadrature);
%             break
        end        
        points_power=points_power+1;
    end
    if points_power>5
        fprintf('Took %d iterations\n',points_power-1)
    end

%     %% Other methods
%     prob=qsimvnauto( sigma_matrix, reachAvoidTubeLB-mu_vector,  reachAvoidTubeUB-mu_vector,myeps,1e5);
%     prob=mvncdf(reachAvoidTubeLB,reachAvoidTubeUB,mu_vector,sigma_matrix);
%     prob=qscmvnv( 50000, sigma_matrix, reachAvoidTubeLB-mu_vector, eye(size(mu_vector,1),size(mu_vector,1)), reachAvoidTubeUB-mu_vector);
end