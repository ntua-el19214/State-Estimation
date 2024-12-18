close all

X_min=-35;
X_max=35;
N=1001;

xValues = linspace(X_min,X_max,N);
pPrior=zeros(N,1);
x0mean=0;
x0std =1;
P_k = x0std;
dx = xValues(2) - xValues(1);

% X_0 distribution
pPrior = x0mean;

X_estimates_EKF = zeros(1,15);
for k=1:15
    y_meas=y(k);
    x_actual = x(k);
    % Prediction step
    x_bar = 0.5*pPrior + bet*pPrior/(1+pPrior.^2) + 8*cos(1.2*k);
      
    P_k_p1 = (bet./(1+pPrior.^2) - 2*bet.*pPrior.^2/(1+pPrior.^2)).^2*P_k + w_std^2;

    % Correction Step
    y_e = y_meas - alpha.*x_bar + x_bar.^2./20;
    S_k_p1 = ((alpha) +x_bar./10)^2.*P_k_p1 + v_std^2;
    K = P_k_p1.*((alpha) +x_bar./10)./S_k_p1;
    x_bar_new = x_bar + K.*y_e;
    X_estimates_EKF(k) = x_bar_new;
    maxAxisScale = 1;
    figure
    quiver(pPrior,0,0,maxAxisScale,  'linewidth',1.5)
    hold on
    quiver(x_bar_new,0,0,maxAxisScale,  'linewidth',1.5)
    quiver(y_meas,0,0,maxAxisScale, 'linewidth',1.5)
    quiver(x_actual,0,0,maxAxisScale, 'linewidth',1.5)
    legend("Prior", "Posterior", "y Measured", "x Actual")
    grid on
    ylabel('Probability Density')
    xlabel('Value')
    titleString = sprintf('x_{%d} Probability Density for \\alpha = %d, \\beta = %d', k, alpha, bet);
    title(titleString);

    %% Save figures
    saveDir = './Figures/EKF'; % Directory to save figures
    if ~exist(saveDir, 'dir') % Check if directory exists
        mkdir(saveDir); % Create directory if it doesn't exist
    end
    filename = sprintf('%s/EKFx_%d_alpha%d_beta%d.png', saveDir, k, alpha, bet);
    % Save the figure as PNG
    saveas(gcf, filename);

    pPrior = x_bar_new;
    P_k_p1 = (1 - K*((alpha) +x_bar./10))*P_k_p1;
end




