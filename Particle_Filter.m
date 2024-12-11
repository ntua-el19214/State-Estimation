rng('default')
%% Parameters
N = 1000;               % Number of particles
num_steps = 15;         % Number of time steps
x_true = x; % True states
y_meas = y; % Measurements

%% Initial State and Particles
x0 = 0;                 % Initial true state
particles = normrnd(x0, 1, [1, N]); % Initialize particles with Gaussian distribution around x0
weights = ones(1, N) / N;           % Initial weights (uniform)

rng('default');
%% Particle Filter Implementation
x_est = zeros(1, num_steps); % Estimated state
for k = 1:num_steps
    % Prediction Step
    for i = 1:N
        w_k = normrnd(0, w_std^2); % Process noise for each particle
        particles(i) = 0.5 * particles(i) + bet * (particles(i) / (1 + particles(i)^2)) + 8 * cos(1.2 * k) + w_k;
    end
    
    % Update Step (Compute Weights)
    for i = 1:N
        % Measurement likelihood (Gaussian likelihood)
        weights(i) = exp(-0.5 * ((y_meas(k) - (alpha * particles(i) + particles(i)^2 / 20)) / v_std)^2);
    end
    
    % Normalize Weights katano
    weights = weights / sum(weights);
    
    % Estimate State (Weighted Average)
    x_est(k) = particles * weights';
    
    % Resampling Step (Systematic Resampling)
    indices = resample(weights, N);
    particles = particles(indices);
    weights = ones(1, N) / N; % Reset weights after resampling
    
    % Plotting the Particle Distribution
    figure;
    hold on
    plot(x_true(k), 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % True state
    plot(x_est(k), 0, 'go', 'MarkerSize', 10, 'LineWidth', 2);  % Estimated state
    plot(linspace(min(x_true), max(x_true), N),particles.*weights, LineWidth=1.5);
    legend('Particles', 'True State', 'Estimated State');
    title(sprintf('Time Step %d', k));
    xlabel('State Value');
    ylabel('Probability');
end

%% Plot True State vs Estimated State
figure;
plot(1:num_steps, x_true, '-r', 'LineWidth', 1.5);
hold on;
plot(1:num_steps, x_est, '-b', 'LineWidth', 1.5);
legend('True State', 'Estimated State');
xlabel('Time Step');
ylabel('State Value');
grid on
title('True State vs Estimated State w/o resampling, N = 50');

%% Resampling algorithm

function indices = resample(weights, N)
    % resampling algorithm
    positions = (0:N-1) / N;
    indices = zeros(1, N);
    cumulative_sum = cumsum(weights);
    i = 1;
    % check if sum of weights at point k is more than the linearly
    % increasing normalized threshold. if not, go to next particle. if yes
    % resample that particle.
    for j = 1:N
        while positions(j) > cumulative_sum(i)
            i = i + 1;
        end
        indices(j) = i;
    end
end
