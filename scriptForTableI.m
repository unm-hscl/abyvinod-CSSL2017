clear
close all
disp('Generate Table 1 (Choose the following options)');
code_flow_flag = input('1-Compute via the scripts for FTBU and DP | 2-Load the matfile\n');
if code_flow_flag == 2
    load('table1_contents.mat');
elseif code_flow_flag == 1
    table_1_generation = 1;
    run_DP_script = 0;
    %%%%% Implementation of FTBU %%%%%%%%%%%
    usePatternSearch=1;
    useFmincon=1;
    displayString_fmincon = 'iter';
    displayString_patternsearch = 'iter';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
else
    error('Invalid option');
end


disp('Comparison of the value functions W_0 (fmincon | patternsearch)');
disp([mean(RA_OP_values,2),mean(RA_OP_PS_values,2)]);
disp('Comparison of the runtimes (fmincon | patternsearch)');
disp([mean(times_FM,2),mean(times_PS,2)]);
