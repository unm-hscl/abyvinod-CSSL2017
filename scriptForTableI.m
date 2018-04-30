clear
close all

table_1_generation = 1;
run_DP_script = 0;
%%%%% Implementation of FTBU %%%%%%%%%%%
usePatternSearch=1;
useFmincon=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scriptForChainDI;

no_of_times=3;
times_FM=[];
times_PS=[];
RA_OP_values=[];
RA_OP_PS_values=[];
for i=1:no_of_times
    scriptForChainDI
    RA_OP_values=[RA_OP_values,RA_OP{40}];
    RA_OP_PS_values=[RA_OP_PS_values,RA_OP_PS{40}];
    times_FM=[times_FM,point_timer_FM];
    times_PS=[times_PS,point_timer_PS];
end

% W_0 value
[mean(RA_OP_values,2),mean(RA_OP_PS_values,2)]
% Runtime
[mean(times_FM,2),mean(times_PS,2)]
