SNR = 0:1:45; % Signal to Noise Ratio (SNR)
N = 1e6; % # of symbols
A = 1; % Maximum aplitude of the cos wave being sent


[bpsk_ser_sim, bpsk_ser_analytic] = BPSK_Simulation(SNR, A, N);
[qpsk_ber_sim, qpsk_ber_analytic, qpsk_ser_sim, qpsk_ser_analytic] = QPSK_Simulation(SNR, A, N);

semilogy(SNR, bpsk_ser_sim, '--x', SNR, bpsk_ser_analytic, '--o', SNR, qpsk_ber_sim, '--', SNR, qpsk_ber_analytic, '--*', SNR, qpsk_ser_sim, '--', SNR, qpsk_ser_analytic, '--*');
ylim([10^-6,10^0])
xlim([0,45])
title('BER & SER: BPSK & QPSK, Theoreticial and Simulation Results over Rayleigh Channel')
legend('M = 2 SER (Simulation)', 'M = 2 SER (Analytical)', 'M = 4 BER (Simulation)', 'M = 4 BER (Analytical)', 'M = 4 SER (Simulation)', 'M = 4 SER (Analytical)');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on


function [ser_simulation, ser_anayltic] = BPSK_Simulation(SNR, A, N)
    absolute_amplitude = A; % Maximum amplitude of the signal
    M = 2; % BPSK 

    % Energy definitions
    Eb = (absolute_amplitude^2) / 2;
    Es = Eb;

    % Voltage of the sent signals
    voltages = zeros(1,M);
    for i = 1:M
        theta = (i-1) *  180;
      
        voltages(i) = -Es * cosd(theta);
    end
      
    % Generate random symbols
    tx = randsrc(1, N, voltages);

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Analytical and simulated error rate arrays
    ser_simulation = zeros(1, length(SNR));
    ser_anayltic = zeros(1, length(SNR));
  

    % Decision alpha
    alpha_opt = 0;

    for i = 1:length(SNR)
        % Rayleigh Fading Channel
        h = (1 / sqrt(2)) * (randn(1, N) + 1i * randn(1, N));

        % Power Calculations
        energy_of_received_symbol = mean(abs(h.^2)) * Es;
        N0 = energy_of_received_symbol / linear_SNR(i);

        % Received Signal with Noise
        rx = h .* tx + sqrt(N0 / 2) * (randn(1, N) + 1i * randn(1, N));

        % Relay Gain for Unit Average Power 
        rx = rx ./ h;

        % Decision Circuit
        decision = real(rx) > alpha_opt;
        decision = voltages(1) * ~decision + voltages(2) * decision;

        % Calculate SER
        ser_simulation(i) = sum(tx ~= decision) / N;% SER Simulated
        ser_anayltic(i) = (1 - (linear_SNR(i) / (1 + linear_SNR(i)))) / 2; % SER Analytical
    end
end

function [ber_simulation, ber_analytical, ser_simulation, ser_analytical] = QPSK_Simulation(SNR, A, N)
    absolute_amplitude = A; % Maximum amplitude of the signal
    M = 4; % QPSK 

    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb * 2;

    % Voltage of the sent signals
    V1 = zeros(1,M);
    V2 = zeros(1,M);

    for i = 1:M
        theta = (180 / M) + (i-1) * 90;
      
        V1(i) = Es * cosd(theta);
        V2(i) = -Es * sind(theta);
    end
    % Generate random symbols
    tx1 = randsrc(1, N, V1);
    tx2 = randsrc(1, N, V2);

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Analytical and simulated error rate arrays
    ber_simulation = zeros(1, length(SNR));
    ber_analytical = zeros(1, length(SNR));
    
    ser_simulation = zeros(1, length(SNR));
    ser_analytical = zeros(1, length(SNR));

    % Decision alpha
    alpha_opt = 0;

    for i = 1:length(SNR)
        % Rayleigh Fading Channel
        h = (1 / sqrt(2)) * (randn(1, N) + 1i * randn(1, N));

        % Power Calculations
        energy_of_received_symbol = mean(abs(h.^2)) * Es;
        N0 = energy_of_received_symbol / linear_SNR(i);

        % Received Signal with Noise
        rx1 = h .* tx1 + sqrt(N0 / 2) * (randn(1, N) + 1i * randn(1, N));
        rx2 = h .* tx2 + sqrt(N0 / 2) * (randn(1, N) + 1i * randn(1, N));

        % Relay Gain for Unit Average Power 
        rx1 = rx1 ./ h;
        rx2 = rx2 ./ h;

        % Decision Circuit
        decision1 = real(rx1) < alpha_opt;
        decision1 = V1(1) * ~decision1 + V1(2) * decision1;
        decision2 = real(rx2) > alpha_opt;
        decision2 = V2(2) * ~decision2 + V2(3) * decision2;

        error1_indices = (tx1 ~= decision1);
        error2_indices = (tx2 ~= decision2);

        % Calculate BER
        ber_simulation(i) = sum(error2_indices) / N; % BER Simulation
        ber_analytical(i) = (1 - (linear_SNR(i) / (1 + linear_SNR(i)))) / 2; % BER Analytical
    
        % Calculate SER
        ser_simulation(i) = (sum(error1_indices | error2_indices)) / N; % SER Simulation
        ser_analytical(i) = 2 * ber_analytical(i) - ber_analytical(i)^2; % SER Analytical
    end
end