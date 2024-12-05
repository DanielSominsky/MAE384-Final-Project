% MAE384 Project %
% Part 4%

% Parameters %

beta0 = 0.3;
A = 5;
gamma = 0.1;
S0 = 990;
I0 = 10;
R0 = 0;
h = 0.1;
T = 30;
t = 0:h:T;

% The Angular Frequencies %

omega1 = 2*pi*365/365;
omega2 = 2*pi*100/365;

%

[S1, I1, R1] = simulate_SIR_variable_beta(t,h,S0,I0,R0,beta0,A,omega1,gamma);
[S2, I2, R2] = simulate_SIR_variable_beta(t,h,S0,I0,R0,beta0,A,omega2,gamma);

[f1, I1_fft] = compute_FFT(I1, h);
[f2, I2_fft] = compute_FFT(I2, h);

% Plots for the time domains and FFT %

figure;

subplot(2, 2, 1);
plot(t, S1,'b', t, I1, 'r', t, R1, 'g');
title('Periodicity: 1 Day');
xlabel('time (days)');
ylabel('population');
legend('S(t)', 'I(t)', 'R(t)');
grid on;

subplot(2, 2, 2);
plot(t, S2, 'b', t, I2, 'r', t, R2, 'g');
title('Periodicity: ~3 Days');
xlabel('time (days)');
ylabel('population');
legend('S(t)', 'I(t)', 'R(t)');
grid on;

subplot(2, 2, 3);
plot(f1, abs(I1_fft));
title('FFT: 1 Day');
xlabel('frequency (Hz)');
ylabel('|FFT| (magnitude)');
grid on;

subplot(2, 2, 4);
plot(f2, abs(I2_fft));
title('FFT: ~3 Days');
xlabel('frequency (Hz)');
ylabel('|FFT| (magnitude)');
grid on;

% Functions %

function [S, I, R] = simulate_SIR_variable_beta(t,h,S0,I0,R0,beta0,A,omega,gamma)

N = length(t);
S = zeros(1, N);
I = zeros(1, N);
R = zeros(1, N);
S(1) = S0;
I(1) = I0;
R(1) = R0;

for i = 1:(N-1)
    beta_t = beta0*(1 + A*sin(omega*t(i)));

    S_i = S(i);
    I_i = I(i);
    R_i = R(i);

    fS = @(S, I) -beta_t .* S .* I/(S0 + I0 + R0);
    fI = @(S, I) beta_t .* S .* I/(S0 + I0 + R0) - gamma .* I;
    fR = @(I) gamma .* I;

    k1_S = h*fS(S_i, I_i);
    k1_I = h*fI(S_i, I_i);
    k1_R = h*fR(I_i);

    k2_S = h*fS(S_i + k1_S/2, I_i + k1_I/2);
    k2_I = h*fI(S_i + k1_S/2, I_i + k1_I/2);
    k2_R = h*fR(I_i + k1_I/2);

    k3_S = h*fS(S_i + k2_S/2, I_i + k2_I/2);
    k3_I = h*fI(S_i + k2_S/2, I_i + k2_I/2);
    k3_R = h*fR(I_i + k2_I/2);

    k4_S = h*fS(S_i + k3_S, I_i + k3_I);
    k4_I = h*fI(S_i + k3_S, I_i + k3_I);
    k4_R = h*fR(I_i + k3_I);

    S(i+1) = S_i + (k1_S + 2*k2_S + 2*k3_S + k4_S) / 6;
    I(i+1) = I_i + (k1_I + 2*k2_I + 2*k3_I + k4_I) / 6;
    R(i+1) = R_i + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6;
end
end

function [frequencies, I_fft] = compute_FFT(I, h)

N = length(I);
I_fft = fft(I);

frequencies = (0:floor(N/2)-1) / (N*h);
I_fft = I_fft(1:floor(N/2));
end

