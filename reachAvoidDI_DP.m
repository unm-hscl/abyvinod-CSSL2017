tic
% Input gridding
udisc=umin:uinc:umax;
no_of_input_points=length(udisc);

% Disturbance gridding
wdisc1D=distmin:distinc:distmax;
DeltaDistArea=distinc*distinc;
wdisc=allcomb(wdisc1D,wdisc1D)';
no_of_dist_points=size(wdisc,2);
pdisc=mvnpdf(wdisc',mean_vector',sigma_matrix)';

% Dynamic programming
Vtplus1=zeros(length(x1),length(x2));
Vtplus1(target_set_LB_indx(1):target_set_UB_indx(1),target_set_LB_indx(2):target_set_UB_indx(2))=ones(target_set_UB_indx(1)-target_set_LB_indx(1)+1,target_set_UB_indx(2)-target_set_LB_indx(2)+1);
Vtplus1=reshape(Vtplus1,[],1);        % change here for dimension
for t=last_time_step-1:-1:0
    Vt=zeros(length(x1)*length(x2),1);
    fprintf('V_t at t=%d...xxxx',t);
    for i=1:length(x1)
        for j=1:length(x2)
            fprintf('\b\b\b\b%1.2f',((i-1)*length(x2)+j)/length(x1)/length(x2));
            xt=[x1(i);x2(j)];
            xtp_no_input=repmat(system_matrix*xt,1,no_of_dist_points)+wdisc;
            EVt_input_vec=zeros(no_of_input_points,1);
%             xtp1=zeros(2,no_of_dist_points);
            for k = 1:no_of_input_points
                xtp1 = xtp_no_input + repmat(input_matrix*udisc(k),1,no_of_dist_points);
                % +1 in the end to get the matlab representation with 1s
                indx_vec=(xtp1-repmat([x1min;x2min],1,no_of_dist_points))./xinc+ones(2,no_of_dist_points);
                % This ensures no spillovers from the boundary
                indx_zeroed=[find(indx_vec(1,:)<1),find(indx_vec(1,:)>length(x1)),find(indx_vec(2,:)<1),find(indx_vec(2,:)>length(x2))];
                % Convert these into indices
                indx_vec(:,indx_zeroed)=NaN;
                indx_vec=round(indx_vec);
                % indx_vec converted to columnwise position
                indx_columnwise=[1 length(x2)]*(indx_vec-[zeros(1,no_of_dist_points);ones(1,no_of_dist_points)]);
                indx_zeroed=find(isnan(indx_columnwise));
                indx_columnwise(indx_zeroed)=1;
                indx=indx_columnwise;
                Vt_temp=Vtplus1(indx,1);
                Vt_temp(indx_zeroed,1)=0;
                EVt_input_vec(k)=pdisc*Vt_temp*DeltaDistArea;
            end
            Vt(i+length(x2)*(j-1),1)=max(EVt_input_vec);
        end        
    end
    fprintf('\b\b\b\bDone!\n');
    Vtplus1=Vt;
end
Vtsquare=reshape(Vt,length(x1),[]);
elapsed_time_DP=toc;
fprintf('Elapsed Time: %1.2f\n',elapsed_time_DP)
