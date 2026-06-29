function [Z1, Y0_hat] = detrend_function(Z0,names,tVec,period_sine,plot_detrend,trend)

if nargin <= 4
    plot_detrend = 1;
    trend = 0;
end    

x = 1:length(Z0);
Period = 4*period_sine;
t = (2*pi/Period).*x;
cycle_sin = sin(t);
cycle_cos = cos(t);

% detrending 
Y0 = Z0;
if period_sine==0
    X0 = [ones(length(Y0),1) (1:length(Z0))'];   % linear trend only
    % X0 = [ones(length(Y0),1) (1:length(Z0))' (1:length(Z0))'.^2];   % linear trend only
else
    X0 = [ones(length(Y0),1) (1:length(Z0))' cycle_sin' cycle_cos'];  % sin and cos
end
if trend==1
    X0 = [ones(length(Y0),1)]; %only intercept
elseif trend ==2
    X0 = [ones(length(Y0),1) (1:length(Z0))' (1:length(Z0))'.^2];
end
betaHat0=(X0'*X0)\(X0'*Y0);
resids0=Y0-X0*betaHat0;
Z1=resids0;
Y0_hat = X0*betaHat0;

Z2 = Z1./std(Z1); %standartized for plotting

k = size(Z2,2);

if plot_detrend ==1
% plot variables and trends
    figure(100)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
    for i = 1:k
        if k<=4
            subplot(2,2,i)
        elseif k==5 || k==6
            subplot(2,3,i)    
        elseif k==7 || k==8
            subplot(2,4,i)
        elseif k==9    
            subplot(3,3,i)
        else 
            subplot(3,4,i)
        end    
        plot(tVec,Z0(:,i),'LineWidth',2,'Color',"r")
        hold on
        plot(tVec,Y0_hat(:,i),'LineWidth',2,'Color',"b",'LineStyle','-.')
        hold off 
        title(names(i))
        set(gca,'fontsize',15)
        grid on
        xlim([tVec(1) tVec(end)])
    end
    % saveas(gcf, 'Figures_money\VariablesAndTrend.jpg')
    
    % plot standard deviations of all variables
    figure(101)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
    for i = 1:k
        if k<=4
            subplot(2,2,i)
        elseif k==5 || k==6
            subplot(2,3,i)   
        elseif k==7 || k==8
            subplot(2,4,i)
        elseif k==9    
            subplot(3,3,i)
        else 
            subplot(3,4,i)
        end    
        plot(tVec,Z2(:,i),'LineWidth',2,'Color',"r")
        title(names(i))
        set(gca,'fontsize',15)
        grid on
        xlim([tVec(1) tVec(end)])
        ylim([-4 4])
        ylabel("Standard deviations")
    end
end    
% saveas(gcf, 'Figures_money\detrendedVariables.jpg')