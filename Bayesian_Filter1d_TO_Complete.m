X_min=-35;
X_max=35;
N=1001;

xValues = linspace(X_min,X_max,N);
pPrior=zeros(N,1);
x0mean=0;
x0std =1;
dx = xValues(2) - xValues(1);

for i=1:N
    pPrior(i) = (1/(2*pi*x0std^2))^0.5  *exp(-(xValues(i)-x0mean)^2/(2*x0std^2));
end


for k=1:10

    y_meas=y(k);
    %correction step
    for x_val_ind=1:N
        x_val_i = xValues(x_val_ind); 

        measurement_model = alpha * x_val_i + (x_val_i^2) / 20;
        likelihood_unnormalized = exp(-((y_meas - measurement_model)^2) / (2 * v_std^2));
        p_posterior_unnormalized (x_val_ind) = pPrior(x_val_ind) * likelihood_unnormalized;
    end
    
    p_posterior = p_posterior_unnormalized/sum(p_posterior_unnormalized*dx);

   figure
   plot(xValues,pPrior)
   hold on
   plot(xValues,p_posterior)
   quiver(y_meas,0,0,0.5)
   quiver(x(k),0,0,0.55)
    
    
    % Prediction Step
    for x_val_new_ind=1:N
        x_val_i_new=xValues(x_val_new_ind);
        for x_val_cur_ind =1:N
            x_val_i_cur = xValues(x_val_cur_ind);
            
            measurement_model_new = alpha * x_val_i_cur + (x_val_i^2) / 20;
            p_Transition_x_to_x_pl  =   (1 / sqrt(2 * pi * w_std^2)) *exp(-1/2*(x_val_i_cur - measurement_model_new)^2/w_std^2);
            p_x_new_to_integrate ( x_val_cur_ind) = p_Transition_x_to_x_pl* p_posterior(x_val_cur_ind);
        end
       
        p_x_new(x_val_new_ind)  = sum(p_x_new_to_integrate)*(xValues(2)-xValues(1));
            
    end
    pPrior=p_x_new;
end 




