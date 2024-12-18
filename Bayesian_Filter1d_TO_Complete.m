close all

X_min=-35;
X_max=35;
N=1001;

xValues = linspace(X_min,X_max,N);
pPrior=zeros(N,1);
x0mean=0;
x0std =1;
dx = xValues(2) - xValues(1);
X_estimates_bayesian = zeros(1,15);
for i=1:N
    pPrior(i) = (1/(2*pi*x0std^2))^0.5  *exp(-(xValues(i)-x0mean)^2/(2*x0std^2));
end


for k=1:15

    y_meas=y(k);
    %correction step
    for x_val_ind=1:N
        x_val_i = xValues(x_val_ind);

        measurement_model = alpha * x_val_i + (x_val_i^2) / 20;
        likelihood_unnormalized = exp(-((y_meas - measurement_model)^2) / (2 * v_std^2));
        p_posterior_unnormalized (x_val_ind) = pPrior(x_val_ind) * likelihood_unnormalized;
    end

    p_posterior = p_posterior_unnormalized/sum(p_posterior_unnormalized)/dx;
    X_estimates_bayesian(k) = (p_posterior*dx*xValues');
    maxAxisScale = max(max(pPrior), max(p_posterior));
    figure
    plot(xValues,pPrior, 'LineWidth',1.5)
    hold on
    plot(xValues,p_posterior, 'LineWidth',1.5)
    quiver(y_meas,0,0,maxAxisScale, 'linewidth',1.5)
    quiver(x(k),0,0,maxAxisScale, 'linewidth',1.5)
    legend("Prior", "Posterior", "y Measured", "x Actual")
    grid on
    ylabel('Probability Density')
    xlabel('Value')
    titleString = sprintf('x_{%d} Probability Density for \\alpha = %d, \\beta = %d', k, alpha, bet);
    title(titleString);

    %% Save figures
    saveDir = './Figures'; % Directory to save figures
    if ~exist(saveDir, 'dir') % Check if directory exists
        mkdir(saveDir); % Create directory if it doesn't exist
    end
    filename = sprintf('%s/x_%d_alpha%d_beta%d.png', saveDir, k, alpha, bet);
    % Save the figure as PNG
    saveas(gcf, filename);

    % Prediction Step
    for x_val_new_ind=1:N
        x_val_i_new=xValues(x_val_new_ind);
        for x_val_cur_ind =1:N
            x_val_i_cur = xValues(x_val_cur_ind);

            measurement_model_new = 0.5*x_val_i_cur+bet*x_val_i_cur/(1+x_val_i_cur^2)+8*cos(1.2*k);
            p_Transition_x_to_x_pl  =   (1 / sqrt(2 * pi * w_std^2)) *exp(-1/2*(x_val_i_new - measurement_model_new)^2/w_std^2);
            p_x_new_to_integrate ( x_val_cur_ind) = p_Transition_x_to_x_pl* p_posterior(x_val_cur_ind);
        end

        p_x_new(x_val_new_ind)  = sum(p_x_new_to_integrate)*dx;
    end
    pPrior=p_x_new;
end




