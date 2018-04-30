tic
    
% Disturbance and state gridding for DP (and FT if underapprox of the value
% function is required)
disp('Creating the grid points for state and disturbance');
wdisc1D=distmin:distinc:distmax;
DeltaDistArea=distinc^state_dimension;
no_of_dist_points=length(wdisc1D)^state_dimension;
no_of_state_points=length(x)^state_dimension;
VtplusDim=[];
reshaping_scaling_matrix=[];
str_iter_vector='allcomb(1:length(x)';
str_wdisc_vector='allcomb(wdisc1D';
for i=1:state_dimension
    if i==1
        reshaping_scaling_matrix=[1];        
    else
        reshaping_scaling_matrix=[reshaping_scaling_matrix,length(x)^(i-1)];
        str_iter_vector=strcat(str_iter_vector,',1:length(x)');
        str_wdisc_vector=strcat(str_wdisc_vector,',wdisc1D');
    end    
end
str_iter_vector=strcat(str_iter_vector,')');
str_wdisc_vector=strcat(str_wdisc_vector,')');
iter_vector=eval(str_iter_vector);
wdisc=eval(str_wdisc_vector);
pdisc=mvnpdf(wdisc,mean_vector',sigma_matrix)';

% Dynamic programming
disp('Entering dynamic programming');
Vtplus1=zeros(no_of_state_points,1);
for indx_iter=1:no_of_state_points
    % Symmetric target
    if all(iter_vector(indx_iter,:)>=target_set_LB_indx(1)) && all(iter_vector(indx_iter,:)<=target_set_UB_indx(1))
        Vtplus1(indx_iter,1)=1;
    end        
end
for t=last_time_step-1:-1:0
    Vt=zeros(no_of_state_points,1);
    fprintf('V_t at t=%d...xxxx',t);
    for indx_iter=1:no_of_state_points
        fprintf('\b\b\b\b%1.2f',indx_iter/no_of_state_points);
% Estimation of time it takes for n=4
%         if indx_iter==1000
%             elapsed_time_DP4=toc;
%             input('Time instant to reach 1000 saved')
%         else
%             fprintf('\n%d',indx_iter);        
%         end
        xt=x(iter_vector(indx_iter,:))';
        xtp_no_input=repmat(system_matrix*xt,1,no_of_dist_points)+wdisc';
        EVt_input_vec=zeros(no_of_input_points,1);
        for k = 1:no_of_input_points
            xtp1 = xtp_no_input + repmat(input_matrix*udisc(k),1,no_of_dist_points);
            % +1 in the end to get the matlab representation with 1s
            indx_vec=(xtp1-repmat(xmin,state_dimension,no_of_dist_points))./xinc+ones(state_dimension,no_of_dist_points);
            % This ensures no spillovers from the boundary
            indx_zeroed=[];
            for j=1:state_dimension
                indx_zeroed=[indx_zeroed,find(indx_vec(j,:)<1),find(indx_vec(j,:)>length(x))];
            end                
            % Convert these into indices
            indx_vec(:,indx_zeroed)=-no_of_state_points;
            indx_vec=round(indx_vec);
            % indx_vec converted to columnwise position
            indx_columnwise=reshaping_scaling_matrix*(indx_vec-[zeros(1,no_of_dist_points);ones(state_dimension-1,no_of_dist_points)]);
            % indx_columnwise is negative iff
            % (find(indx_vec(j,:)<1),find(indx_vec(j,:)>length(x))) or
            % indx
            indx_zeroed=find(indx_columnwise<0);
            indx_columnwise(indx_zeroed)=1;
            indx=indx_columnwise;
            Vt_temp=Vtplus1(indx,1);
            Vt_temp(indx_zeroed,1)=0;
            EVt_input_vec(k)=pdisc*Vt_temp*DeltaDistArea;
        end
        Vt(indx_iter,1)=max(EVt_input_vec);
    end
    fprintf('\b\b\b\bDone!\n');
    Vtplus1=Vt;
end
% Vtsquare=reshape(Vt,length(x),[]);
eval(sprintf('terminal_RA_prob{%d}=Vt;',state_dimension));
testing_elapsed_time_DP=toc;
eval(sprintf('elapsed_time_DP{%d}=testing_elapsed_time_DP;',state_dimension));
fprintf('Elapsed Time: %1.2f\n',testing_elapsed_time_DP)
