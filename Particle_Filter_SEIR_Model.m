rng('default')
close all
numOfSteps = 11;
k_start = 100;
%% Initialize particles & weights
N = 3000;
particles_s = zeros(1, N);
particles_e = zeros(1, N);
particles_i = zeros(1, N);
% setting covariance of s and e to be 0.8 as the two distributions are
% added
for iParticle = 1:N
    particles_s(i) = S(k_start) -1*rand(1)*0.01;
    particles_e(i) = E(k_start) + 0 ;
    particles_i(i) = I(k_start)+rand(1)*0.001;
end

weights = ones(1, N) / N;           % Initial weights (uniform)

S_est = ones(1, numOfSteps).*S(k_start) -1*rand(1)*0.01; 
E_est = ones(1, numOfSteps).*E(k_start) +rand(1)*0.01;
I_est = ones(1, numOfSteps).*I(k_start)+rand(1)*0.01;

%% Particle filter implementation
for k = 1:numOfSteps
    % Prediction step for particle filter
    y_k = y(k) ;
    % Update step: calculate weights based on observation likelihood
    for i = 1:N
        predicted_measurement = 0.2 * particles_e(i) + particles_i(i);
        weights(i) = (predicted_measurement / y_k) * exp(-0.5 * ((log(y_k) - log(predicted_measurement)) / 0.3)^2);
    end

    % Normalize weights
    weights = weights / sum(weights);

    % Estimate state
    S_est(k) = particles_s * weights';
    E_est(k) = particles_e * weights';
    I_est(k) = particles_i * weights';

    % Resample particles
    indices = resample(weights, N);
    particles_s = particles_s(indices);
    particles_e = particles_e(indices);
    particles_i = particles_i(indices);
    
    for i = 1:N
        w1_k = rand(1)*0.5;
        w2_k = rand(1)*0.3;
        w3_k = rand(1)*0.3;
        
        particles = SEIR_Dynamics(particles_s(i), particles_e(i), particles_i(i), w1_k, w2_k, w3_k);

        particles_s(i) = particles(1);
        particles_e(i) = particles(2);
        particles_i(i) = particles(3);
    end 
    % Reset weights after resampling
    weights = ones(1, N) / N;
end

%% Plot results for E
figure;
plot(k_start:k_start+numOfSteps-1, E(k_start:k_start+numOfSteps-1), '-r', 'LineWidth', 1.5);
hold on;
plot(k_start:k_start+numOfSteps-1, E_est, '-b', 'LineWidth', 1.5);
legend('True E State', 'Estimated E State');
xlabel('Time Step');
ylabel('State Value');
% ylim([0 0.01])
grid on;
title('True E State vs Estimated E State with Particle Filter');

%% Plot results for S
figure;
plot(k_start:k_start+numOfSteps-1, S(k_start:k_start+numOfSteps-1), '-r', 'LineWidth', 1.5);
hold on;
plot(k_start:k_start+numOfSteps-1, S_est, '-b', 'LineWidth', 1.5);
legend('True S State', 'Estimated S State');
xlabel('Time Step');
ylabel('State Value');
% ylim([0.9 1])
grid on;
title('True S State vs Estimated S State with Particle Filter');

%% Plot results for I
figure;
plot(k_start:k_start+numOfSteps-1, I(k_start:k_start+numOfSteps-1), '-r', 'LineWidth', 1.5);
hold on;
plot(k_start:k_start+numOfSteps-1, I_est, '-b', 'LineWidth', 1.5);
legend('True I State', 'Estimated I State');
xlabel('Time Step');
ylabel('State Value');
% ylim([0 0.01])
grid on;
title('True I State vs Estimated I State with Particle Filter');

%% Resampling function
function indices = resample(weights, N)
    % Systematic resampling algorithm
    positions = (rand + (0:N-1)) / N;
    indices = zeros(1, N);
    cumulative_sum = cumsum(weights);
    i = 1;
    for j = 1:N
        while positions(j) > cumulative_sum(i)
            i = i + 1;
        end
        indices(j) = i;
    end
end

fprintf("MSE for k = %f\n", k_start)
mse_s = sum((S(k_start:k_start+numOfSteps-1) - S_est).^2)
mse_e = sum((E(k_start:k_start+numOfSteps-1) - E_est).^2)
mse_i = sum((I(k_start:k_start+numOfSteps-1) - I_est).^2)

           %K = 1        K = 50     K = 100  
MSE_est = [0.0018        5.5163e-04  0.0269     ;  % S
           9.7774e-04    4.7689e-04  5.9483e-04 ;  % E
           4.7082e-04    3.3833e-04  1.8591e-04 ]; % I

mean(MSE_est,2)

           % K = 1       K = 50               K = 100
MSE_act = [3.8459e-04    0.0064             0.0011;  % S
           1.2344e-06    0.0012             0.0020;  % E
           4.7865e-06    1.3325e-04         6.3389e-04]; % I
mean(MSE_act,2)