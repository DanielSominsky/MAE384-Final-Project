%Initial Parameters
h = 1;                %Time step (in days)
T = 100;              %Total simulation time (in days)
N = 1000;             %Total population
t = 0:h:T;            %Time vector

%Initial Conditions
S0 = 990;             %Initial susceptible people
I0 = 10;              %Initial infected people
R0 = 0;               %Initial recovered people

%Parameter Combinations: [beta, gamma]
param_sets = [
    0.3, 0.1;         %Seasonal Influenza
    1.0, 0.1;         %COVID-19
    2.0, 0.2          %Measles
];
param_names = {'Seasonal Influenza', 'COVID-19', 'Measles'};

%Solve SIR model for each parameter set
figure;
for scenario = 1:size(param_sets, 1)
    beta = param_sets(scenario, 1);
    gamma = param_sets(scenario, 2);

    %Initialize State Vars
    S = zeros(1, length(t));
    I = zeros(1, length(t));
    R = zeros(1, length(t));
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;

    %Runge-Kutta 4th Order
    for i = 1:(length(t)-1)
        %Current state
        S_i = S(i);
        I_i = I(i);
        R_i = R(i);

        %Define Derivatives
        f1_S = @(S, I) -beta / N * S * I;
        f1_I = @(S, I) beta / N * S * I - gamma * I;
        f1_R = @(I) gamma * I;

        %RK4 steps for S
        k1_S = h * f1_S(S_i, I_i);
        k1_I = h * f1_I(S_i, I_i);
        k1_R = h * f1_R(I_i);

        k2_S = h * f1_S(S_i + k1_S / 2, I_i + k1_I / 2);
        k2_I = h * f1_I(S_i + k1_S / 2, I_i + k1_I / 2);
        k2_R = h * f1_R(I_i + k1_I / 2);

        k3_S = h * f1_S(S_i + k2_S / 2, I_i + k2_I / 2);
        k3_I = h * f1_I(S_i + k2_S / 2, I_i + k2_I / 2);
        k3_R = h * f1_R(I_i + k2_I / 2);

        k4_S = h * f1_S(S_i + k3_S, I_i + k3_I);
        k4_I = h * f1_I(S_i + k3_S, I_i + k3_I);
        k4_R = h * f1_R(I_i + k3_I);

        %Update Values
        S(i+1) = S_i + (k1_S + 2*k2_S + 2*k3_S + k4_S) / 6;
        I(i+1) = I_i + (k1_I + 2*k2_I + 2*k3_I + k4_I) / 6;
        R(i+1) = R_i + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6;
    end

    %Plotting
    subplot(3, 1, scenario);
    plot(t, S, 'b', 'LineWidth', 1.5, 'DisplayName', 'Susceptible');
    hold on;
    plot(t, I, 'r', 'LineWidth', 1.5, 'DisplayName', 'Infected');
    plot(t, R, 'g', 'LineWidth', 1.5, 'DisplayName', 'Recovered');
    title(['SIR Model: ', param_names{scenario}, ' (\beta=', num2str(beta), ', \gamma=', num2str(gamma), ')']);
    xlabel('Time (days)');
    ylabel('Population');
    legend('Location', 'best');
    grid on;
end
