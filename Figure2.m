clear
close all

%%%%% Implementation of FTBU %%%%%%%%%%%
usePatternSearch=1;
useFmincon=1;
displayString_fmincon = 'off';
displayString_patternsearch = 'iter';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generate Figure 2 (Choose the following options)');
code_flow_flag = input('1-Compute via the scripts for FTBU and DP | 2-Load the matfile\n');
if code_flow_flag == 2
    load('Figure2_0x05_onX_bothPSandFM');
    run_DP_script = 0;
elseif code_flow_flag == 1
    scriptForComparison
else
    error('Invalid option');
end
show_title = 1;

%% Comment the next line and run this after running scriptForComparison.m
% This simulation was run on my laptop (which has slightly lower configs
% than the lab desktop) --- hence the higher computational runtimes

figure(2)
clf
colorbar;
hold on;
surf(x1,x2,Vtsquare');
minc=0;
maxc=1;
set(gca,'clim',[minc maxc]);
set(gca,'Fontsize',15);
axis([x1min x1max x2min x2max minc maxc]);
if show_title == 1
    title('ReachAvoid DP approximation');
end
xlabel('$x_1$','interpreter','latex','Fontsize',20);
ylabel('$x_2$','interpreter','latex','Fontsize',20);
box on;
grid on;
view([0,90]);

if useFmincon==1
    figure(3)
    clf
    colorbar;
    hold on;
    minc=0;
    maxc=1;
    % Plotting assumes 10 time steps
    surf(x1,x2,RA_OP');
    set(gca,'clim',[minc maxc]);
    set(gca,'Fontsize',15);
    if show_title == 1
        title('ReachAvoid Underapproximation --- fmincon');
    end
    axis([x1min x1max x2min x2max minc maxc]);
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    view([0,90]);

    %% Relative error computation
    consvPos=consv;
    consvPos(Vtsquare<0.01)=NaN;
    relerror=consvPos./Vtsquare*100;
    figure(5)
    clf
    minc=0;
    maxc=100;
    surf(x1,x2,relerror')
    colorbar;
    set(gca,'clim',[minc maxc]);
    set(gca,'Fontsize',15);
    if show_title == 1
        title('Relative error in conservativeness --- fmincon vs DPBDA');
    end
    axis([x1min x1max x2min x2max minc maxc]);
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    view([0,90]);
end

if usePatternSearch==1
    figure(4)
    clf
    colorbar;
    hold on;
    minc=0;
    maxc=1;
    % Plotting assumes 10 time steps
    surf(x1,x2,RA_OP_PS');
    set(gca,'clim',[minc maxc]);
    set(gca,'Fontsize',15);
    if show_title == 1
        title('ReachAvoid Underapproximation --- patternsearch');
    end
    axis([x1min x1max x2min x2max minc maxc]);
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    view([0,90]);

    %% Relative error computation
    consvPos=consv_PS;
    consvPos(Vtsquare<0.01)=NaN;
    relerror=consvPos./Vtsquare*100;
    figure(6)
    clf
    minc=0;
    maxc=100;
    surf(x1,x2,relerror')
    colorbar;
    set(gca,'clim',[minc maxc]);
    set(gca,'Fontsize',15);
    if show_title == 1
        title('Relative error in conservativeness --- patternsearch vs DPBDA');
    end
    axis([x1min x1max x2min x2max minc maxc]);
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    view([0,90]);

    %% 30% or less relative error
    figure(7)
    clf
    indx_relerror=relerror<30;
    contourf(x1,x2,indx_relerror',[1 1]);
    if show_title == 1
        title('Filled contour for <30% error --- patternsearch');
    end
    legend('Relative error <30%')
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    axis([x1min x1max x2min x2max minc maxc]);
    set(gca,'FontSize',15);
    view([0,90]);
end

if useFmincon==1 && usePatternSearch==1
    diffPSminusFM=consv-consv_PS;
    diffPSminusFM(Vtsquare<0.01)=NaN;
    figure(8);
    clf
    minc=min(min(diffPSminusFM));
    maxc=max(max(diffPSminusFM));
    surf(x1,x2,diffPSminusFM')
    colorbar;
    set(gca,'clim',[minc maxc]);
    set(gca,'Fontsize',15);
    if show_title == 1
        title('Relative error in conservativeness --- patternsearch vs fmincon');
    end
    axis([x1min x1max x2min x2max minc maxc]);
    xlabel('$x_1$','interpreter','latex','Fontsize',20);
    ylabel('$x_2$','interpreter','latex','Fontsize',20);
    box on;
    grid on;
    view([0,90]);
end
