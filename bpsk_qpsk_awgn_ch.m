SNR = 0:1:16; % Signal to Noise Ratio (SNR)
N = 1e6; % # of symbols
A = 1; % Maximum aplitude of the cos wave being sent
format long;

[qpsk_ber_sim, qpsk_ber_analytic, qpsk_ser_sim, qpsk_ser_analytic] = QPSK_Simulation(SNR,A, N);
[bpsk_ser_sim, bpsk_ser_analytic] = BPSK_Simulation(SNR, A, N);

semilogy(SNR, bpsk_ser_sim, '--x', SNR, bpsk_ser_analytic, '--o', SNR, qpsk_ber_sim, '--*', SNR, qpsk_ber_analytic, SNR, qpsk_ser_sim, '--*', SNR, qpsk_ser_analytic);
ylim([10^-6,10^0])
xlim([0,13])
title('BER & SER: BPSK & QPSK, Theoreticial and Simulation Results Over AWGN Channel')
legend('M = 2 SER (Simulation)', 'M = 2 SER (Analytical)', 'M = 4 BER (Simulation)','M = 4 BER (Analytical)', 'M = 4 SER (Simulation)','M = 4 SER (Analytical)');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on

function [ser_simulation, ser_analytical] = BPSK_Simulation(SNR, A, N)
    absolute_amplitude = A; % Absolute value of the amplitudes of the signals being sent
    M = 2; % BPSK

    % Values of the noiseless voltages collected by the decision circuit
    voltages = zeros(1,M); 

    % Energy definitions
    symbol_energy = (absolute_amplitude ^ 2) / 2; 
    Es = symbol_energy; 

    for i = 1:M
        theta = (i-1) *  180;
      
        voltages(i) = symbol_energy * cosd(theta);
    end
    
    % Randomly generating the symbols
    tx = randsrc(1, N, voltages);

    % Power Calculations
    linear_SNR = 10 .^ (SNR./10);     

    N0 = Es ./ linear_SNR;

    signalpower = (symbol_energy*Es) / 2;
    
    % Received signals
    rx = zeros(length(SNR),N);
    for i = 1:length(SNR) 
        rx(i,:) = awgn(tx, linear_SNR(i), signalpower, 'linear');
    end
    
    alpha_opt = 0;

    % Buffer keeping the decisions for every SNR value
    buffer1 = zeros(1,N);
    
    % Array keeping the decisions
    decision1 = zeros(length(SNR), N);

    % Analytical and simulated error rate arrays
    ser_simulation = zeros(1,length(SNR));
    ser_analytical = zeros(1,length(SNR));

    % Decision circuit
    for i = 1:length(SNR)                           
    buffer1(rx(i,:) > alpha_opt) = voltages(1);
    buffer1(rx(i,:) < alpha_opt) = voltages(2);

    decision1(i,:) = buffer1;
    
    ser_simulation(i) = symerr(tx, decision1(i,:)) / N; % SER Simulation
    ser_analytical(i) = qfunc(sqrt(2 * Es / N0(i))); % SER Analytical

    end
end
    

function [ber_simulation, ber_analytical, ser_simulation, ser_analytical] = QPSK_Simulation(SNR, A, N)
    absolute_amplitude = A; % Absolute value of the amplitudes of the signals being sent
    M = 4; % QPSK

    % Values of the noiseless voltages collected by the decision circuit
    V1 = zeros(1,M); 
    V2 = zeros(1,M);

    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb * 2;

    for i = 1:M
        theta = (180 / M) + (i-1) * 90;
      
        V1(i) = Es * cosd(theta);
        V2(i) = -Es * sind(theta);
    end

    % Randomly generating the symbols (tx1 and tx2 parts of the symbol tx) 
    tx1 = randsrc(1, N, V1);
    tx2 = randsrc(1, N, V2);

    % Power Calculations
    linear_SNR = 10 .^ (SNR./10);     

    N0 = Es ./ linear_SNR;

    signalpower = (Eb*Es) / 2;
    
    % Received signals
    rx1 = zeros(length(SNR),N);
    rx2 = zeros(length(SNR),N);
    for i = 1:length(SNR) 
        rx1(i,:) = awgn(tx1, linear_SNR(i), signalpower, 'linear');
        rx2(i,:) = awgn(tx2, linear_SNR(i), signalpower, 'linear');
    end
    
    alpha_opt = 0;

    % Buffers keeping the decisions for every SNR value
    buffer1 = zeros(1,N);
    buffer2 = zeros(1,N);
    
    % Array keeping the decisions
    decision1 = zeros(length(SNR), N);
    decision2 = zeros(length(SNR), N);

    % Analytical and simulated error rate arrays
    ber_simulation = zeros(1,length(SNR));
    ber_analytical = zeros(1, length(SNR));
    ser_simulation = zeros(1, length(SNR));
    ser_analytical = zeros(1,length(SNR));

    % Decision circuit
    for i = 1:length(SNR)                           
        buffer1(rx1(i,:) > alpha_opt) = V1(1);
        buffer1(rx1(i,:) < alpha_opt) = V1(2);
        buffer2(rx2(i,:) > alpha_opt) = V2(3);
        buffer2(rx2(i,:) < alpha_opt) = V2(2);
    
        decision1(i,:) = buffer1;
        decision2(i,:) = buffer2;

        error1_indices = (tx1 ~= decision1(i,:));
        error2_indices = (tx2 ~= decision2(i,:));
    
        ber_simulation(i) = sum(error2_indices) / N; % BER Simulation
        ber_analytical(i) = qfunc(sqrt(2 * Es / N0(i))); % BER Analytical
         
        ser_simulation(i) =  sum(error1_indices | error2_indices) / N; % SER Simulation
        ser_analytical(i) = 2 * qfunc(sqrt(2 * Es/N0(i))) - qfunc(sqrt(2 * Es/N0(i)))^2; % SER Analytical

    end
    
end