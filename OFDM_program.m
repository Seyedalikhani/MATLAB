clc
clear
close all

% (((((((((((((((((((( Main Parametes for Programe ))))))))))))))))))))

ML=0;   % You need to set one parameter to 1 in others must be 0
MAP=0;
Coding=1;

% number of loops (messages to be sent for each carrier)
Nmes = 48; % change to higher number later when one iteration of the loop works!

N0=2000;           % Noise power

figure_flag=0;     % you can change it to 1 if you one to see figures in the output
if figure_flag==1
    Nmes=1;
end


% (((((((((((((((((((( Main Parametes for Programe ))))))))))))))))))))

format long
%% System parameters
K = 64; % Number of subcarriers
N = 64; % number of samples equal to K.
L = floor(0.1*K); % cyclic prefix samples "approx 10%"
Ts = 0.001; % signaling time
Tobs = (N*Ts)/(N+L); % observation time
Tcp = Ts - Tobs; % Length of cyclic prefix
fDelta = 1/Tobs; % separation of carriers
fsamp = N*fDelta; % Choosen
fsampMin = K*fDelta; % The lowest rate we can sample according to theory.
sampleTimes = 0:Ts/(N+L):Ts - Ts/(N+L); % times where the N samples are found (including cyclic prefix).
N_analog = 1000*(N+L); %number of samples per Ts for the simulated analog signal.
sampleIndexes = 1:floor(N_analog/(L+N)):(L+N)*floor(N_analog/(L+N)); % sample times
% indexes for the analog signal where the N samples are.



fc = (100)*N_analog; % Carrier frequency

% Channel model: sum_k channelGains(k) \delta(t - channelTaus(k))

Nchannel = 4; % number of taps or time delays in the channel.
% Preciseply: length(channelTaus).

channelGains = [0.9, 0.3 0.03 0.02];
channelTaus = [0 0.02 0.3 0.33]*Tcp; % We have a very generous choice of Tcp,


% we see that Th is approximately (1/3)Tcp. But let's be safe and we chose
% Tcp approximately 10% of Tobs.

channelSampleTimes = floor(channelTaus*(N_analog/Ts)); % the sample times for
% the analog signal correspoinding to channelTaus.

SamplesDuringTcp = floor(Tcp*(N_analog/Ts)); % simply the amount of samples
% of the analog signal during a time interval of length Tcp.




%% We produce a bunch of subcarrierers or subchannels (K such to be precise)
subcarrier = cell(K,1);

for k = 1:K
    % 64-QAM for all of them
    carrier.k = 6;
    carrier.M = 2^carrier.k; %sqrt(M) must be integer!
    carrier.sqrtM = sqrt(carrier.M);


    % Constellation same for all, Dmin = 2.
    carrier.Ichoices = -carrier.sqrtM+1:2:carrier.sqrtM-1; % x-coordinates
    carrier.Qchoices = -carrier.sqrtM+1:2:carrier.sqrtM-1; % y-coordinates

    % Constellation matrix
    carrier.C = ones(carrier.sqrtM,1)*carrier.Ichoices + carrier.Qchoices'*ones(1,carrier.sqrtM)*1i;

    % Priors
    % carrier.p = rand(1,carrier.M); random prior
    %carrier.p = (1:carrier.M); % some specific choice of prior
    carrier.p = ones(1,carrier.M); % flat prior
    carrier.p = carrier.p/norm(carrier.p,1); % sum over the prior should be 1.
    carrier.P = reshape(carrier.p, [carrier.sqrtM, carrier.sqrtM]); % Prior in matrix form;


    % Cumulative distribution
    carrier.pCum = carrier.p*tril(ones(carrier.M, carrier.M))';

    % Shifted cumulative distribution
    carrier.pCumShift = [0, carrier.pCum(1:carrier.M-1)];

    % save to array of subcarriers
    subcarrier{k,1} = carrier;
end




% Generating Default Subcarrier and Message Matrixes
M=sqrt(K);
SubCarrier_Mat=zeros(M,M);
Message_Mat=zeros(M,M);
for m=1:M
    for n=1:M
        SubCarrier_Mat(n,m)=-M+1+2*(m-1)-(-M+1+2*(n-1))*1i;
        Message_Mat(n,m)=n+8*(m-1);
    end
end




%% We start the main loop


% These guys are increased in the symbol detection later.
NrSymbolErrorsML = 0;
NrSymbolErrorsMAP = 0;
NrSymbolErrorsCoding = 0;

for iterations = 1:Nmes

    %% Transmitter

    % choose messages and symbols
    messages = zeros(K,1);
    a = zeros(K,1);

    % create the random symbol to be sent for each carrier and put it in the vector a
    if ML==1
        for k = 1:K
            % Draw a message from the prior
            M = subcarrier{k,1}.M;
            probVec = rand*ones(1,M);
            messageChoosen = find((ones(1,M) - sign((probVec - ...
                subcarrier{k,1}.pCum).*(probVec - subcarrier{k,1}.pCumShift))));
            messages(k) = messageChoosen;
            a(k) = subcarrier{k,1}.C(messageChoosen);
        end
    end
    if MAP==1
        % MAP Message Generation
        while (sum(a)==0 || length(a)~=64)
            a=MAP_Generator(SubCarrier_Mat)';
        end
    end
    % the a-vector contains the symbols to be sent.


    if Coding==1
        if mod(iterations,3)~=0
            A_Message(iterations).Vec=randperm(K);
            a=SubCarrier_Mat(reshape(A_Message(iterations).Vec,64,1));
            A_Message(iterations).Message=reshape(A_Message(iterations).Vec,8,8);
            A_Message(iterations).a=a;
            A_Message(iterations).Type='Main Message';
        elseif mod(iterations,3)==0
            ind1=(iterations/3-1)*3+1;
            ind2=ind1+1;
            M1=[];
            M2=[];
            redundant_M=[];
            p = randperm(K);
            for m=1:64
                M1=[M1 Message_Mat(find(Message_Mat== A_Message(ind1).Vec(m)))];
                M2=[M2 Message_Mat(find(Message_Mat== A_Message(ind2).Vec(m)))];
                redundant_M=[redundant_M p(mod(M1(m) + M2(m), 64) + 1)];
            end
            A_Message(iterations).Vec=redundant_M;
            A_Message(iterations).Message=reshape(redundant_M,8,8);
            A_Message(iterations).p=p;
            a=SubCarrier_Mat(reshape(A_Message(iterations).Message,64,1));
            A_Message(iterations).a=a;
            A_Message(iterations).Type='Redundancy';
        end
    end



    % Construct the matrix Qt and the matrix Qt'
    krc=(K-2)/2;
    Q(N-krc+1:N,1:krc)=eye(krc);
    Q(1:K-krc,krc+1:K)=eye(K-krc);
    Qt=Q;
    Qt_prime=Q';

    % Vector gk and fk k=0,...,K-1
    if mod(K,2)==0
        frc=fc-fDelta/2;
        g=(0:K-1)-(K-2)/2;
        fk=frc+g*fDelta;
    else
        frc=fc;
        g=(0:K-1)-(K-1)/2;
        fk=frc+g*fDelta;
    end


    % IFFT Samples in Tobs: [Tcp - Ts]  size (N*1)
    X1=ifft(N*Qt*a,N);
    % X1= (x[n]=x(nTobs/N)=sigma k=0,...K-1 (ak*exp(j*2*pi*gk*n*/N)),
    % x[0]=X1[1]=sigma k=0,...K-1 (ak),  X1[1]=sum(a)

    % Add the cyclic prefix
    X=[X1(N-L+1:N); X1];    % [0-Ts]


    t = linspace(0,Ts, N_analog);
    X_Re=real(X);      % Real part of baseband OFDM samples
    X_Im=imag(X);      % Imaginary part of baseband OFDM samples
    M1=13;
    X_Re_signal = zeros(size(t));   % Real part of baseband OFDM analog signal
    X_Im_signal = zeros(size(t));   % Imaginary part of baseband OFDM analog signal
    for n = 1:L+N
        % Compute sinc function values for each sample
        X_Re_signal = X_Re_signal + X_Re(n) * truncated_sinc(t - (n-1) / fsamp, fsamp, M1);
        X_Im_signal = X_Im_signal + X_Im(n) * truncated_sinc(t - (n-1) / fsamp, fsamp, M1);
    end

    if figure_flag==1
        figure(1);
        ax(1)=subplot(2,1,1);
        plot(t, X_Re_signal, 'b', 'LineWidth', 1.5); hold on;
        stem((0:L+N-1) * Tobs/N, X_Re, 'r', 'LineWidth', 1);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Real part of baseband Signal and OFDM samples (X-Re(t))');
        legend('Real Baseband Signal', 'Original Samples');
        grid on;


        ax(2)=subplot(2,1,2);
        plot(t, X_Im_signal, 'b', 'LineWidth', 1.5); hold on;
        stem((0:L+N-1) * Tobs/N, X_Im, 'r', 'LineWidth', 1);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Imaginary part of baseband Signal and OFDM Samples (X-Im(t))');
        legend('Imaginary Baseband Signal', 'Original Samples');
        grid on
    end



    % Modulated Signal: Z(t)=x_I(t)cos(2*pi*fc*t)-x_Q(t)cos(2*pi*fc*t)
    Y_Bandpass=X_Re_signal.*cos(2*pi*fc*t)-X_Im_signal.*sin(2*pi*fc*t);
    if figure_flag==1
        figure(2)
        plot(t,Y_Bandpass)
        title('High frequency OFDM signal (y(t)=X-Re(t)cos(2*pi*fc*t)-X-Im(t)sin(2*pi*fc*t))');
        grid on
    end


    %% Channel and Noise
    % h[n]=sigma n=0,..,3 a_n*delta(t-tha_n)
    % Y_Bandpass*h[n]=sigma n=0,..,3 a_n*Y_Bandpass(t-tha_n)

    % Time Domain BP signal with multipatch
    Y_Bandpass_Multipath=zeros(1,N_analog);

    for j=1:Nchannel
        Taus_index=find(abs(t-channelTaus(j))==min(abs(t-channelTaus(j))));
        ch_ray=[zeros(1,Taus_index-1) Y_Bandpass(1:N_analog-Taus_index+1)];
        ch_ray=channelGains(j)*ch_ray;
        Y_Bandpass_Multipath=Y_Bandpass_Multipath+ch_ray;
        if figure_flag==1
            figure(3);
            ax(j)=subplot(Nchannel,1,j);
            plot(t, ch_ray, 'b', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Amplitude');
            if j==1
                title('4-ray channel multipath channel ')
            end
            grid on
        end
    end


    % Channel Frequency Reponse
    N_H = K;              % Number of frequency points
    % Frequency response calculation
    H = zeros(1, N_H);       % Initialize the frequency response

    for c = 1:Nchannel
        H = H + channelGains(c) * exp(-1j * 2 * pi * fk * channelTaus(c));
    end

    % Magnitude and phase
    magnitude = abs(H);
    phase = angle(H);

    % Frequency response of the channel during the BW of the OFDM signal
    H_f=magnitude.*exp(1i*phase);


    if figure_flag==1
        figure(4)
        % Magnitude response
        subplot(2, 1, 1);
        plot(fk / 1e6, 20 * log10(magnitude),'*'); % Convert frequency to MHz and magnitude to dB
        title('Magnitude Response of Multipath Channel');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude (dB)');
        grid on;

        % Phase response
        subplot(2, 1, 2);
        plot(fk / 1e6, unwrap(phase),'*'); % Unwrapped phase
        title('Phase Response of Multipath Channel');
        xlabel('Frequency (MHz)');
        ylabel('Phase (radians)');
        grid on;


        if figure_flag==1
            figure(5);
            plot(t,Y_Bandpass_Multipath)
            grid on
            xlabel('Time (s)');
            ylabel('Amplitude');
            title('Bandpass signal in receiver assuming multipath channel ');
        end
    end



    % Filtered Modulated Signal and Receiver
    N1 = 1024;
    t1 = -N1/2:N1/2-1; % Time vector for sinc
    fcut=5000;

    % Create the sinc filter
    sinc_filter = sinc(2 * fcut * t1 / fsamp);
    % Apply a window to truncate (e.g., Hamming window)
    window = hamming(length(sinc_filter))';
    truncated_sinc_filter = sinc_filter .* window;
    truncated_filter = truncated_sinc_filter / sum(truncated_sinc_filter); % Normalize filter

    if figure_flag==1
        figure(6)
        subplot(2,1,1)
        plot(t1,truncated_filter)
        title('Sinc low pass filter')
        grid on
    end

    %Compute FFT
    sinc_fft = fftshift(fft(truncated_filter, N1)); % Center zero frequency
    frequencies = (-N1/2:N1/2-1) * (fsamp / N1);    % Frequency axis

    %Normalize FFT magnitude
    sinc_fft_magnitude = abs(sinc_fft) / max(abs(sinc_fft));

    if figure_flag==1
        subplot(2,1,2)
        plot(frequencies, sinc_fft_magnitude);
        title('FFT of Sinc Function');
        xlabel('Frequency (Hz)');
        ylabel('Normalized Magnitude');
    end


    % Convolve the signal with the truncated sinc filter
    Zlp_real = conv(Y_Bandpass_Multipath.*cos(2*pi*fc*t), truncated_filter, 'same');
    Zlp_imaginary = conv(Y_Bandpass_Multipath.*(-sin(2*pi*fc*t)), truncated_filter, 'same');


    if figure_flag==1
        figure(7)
        subplot(2, 1, 1);
        plot(t, X_Re_signal, 'b', 'LineWidth', 1.5); hold on
        plot(t, Zlp_real,'r', 'LineWidth', 1.5);
        title('Real part of baseband signal after filtering')
        legend('Transmitter','Receiver');
        xlabel('Time');
        ylabel('Amplitude');
        grid on;

        subplot(2, 1, 2);
        plot(t, X_Im_signal, 'b', 'LineWidth', 1.5); hold on
        plot(t, Zlp_imaginary, 'r','LineWidth', 1.5);
        title('Imaginary part of baseband signal after filtering')
        legend('Transmitter','Receiver');
        xlabel('Time');
        ylabel('Amplitude');
        grid on
    end


    Zlp=(Zlp_real(sampleIndexes)-Zlp_imaginary(sampleIndexes)*1i);
    Zlp=Zlp(L+1:N+L)';
    Heq=2;
    R=fft(Zlp,N)*Heq;

    % Add Noise
    sigma = sqrt(N0 / 2);
    noise_real = sigma * randn(1, K);
    noise_imag = sigma * randn(1, K);
    R_real=real(R)+noise_real';
    R_imag=imag(R)+noise_imag';
    R=R_real+R_imag*1i;

    a_prime=1/N*Qt_prime*R./conj(H_f');

    Heq_distance=[
        1.8000    7.0945
        1.8500    6.2156
        1.9000    5.7727
        1.9500    5.3383
        2.0000    5.2840
        2.0500    5.4701
        2.1000    6.0299
        2.1500    6.6138
        2.2000    7.5416];
    if figure_flag==1
        figure(8)
        plot(Heq_distance(:,1),Heq_distance(:,2),'--*')
        title('OFDM Average Symbol Distance v.s. Heq in the Reciever, Heq(optimum)=2')
        grid on
    end



    if figure_flag==1
        figure(9)
        for j=-6:2:6
            plot([-8,8],[j,j],LineWidth=2,Color="k")
            hold on
            plot([j,j],[-8,8],LineWidth=2,Color="k")
        end
    end

    Decoded=[];
    l=1:K;
    P=2*l/(K*(K+1));
    for s=1:K
        distance_vector=[];
        Certainly_NoError=0;
        for r=1:K
            if ML==1
                distance_vector=[distance_vector norm(SubCarrier_Mat(r)-a_prime(s))];
            end
            if MAP==1
                distance_vector=[distance_vector norm(SubCarrier_Mat(r)-a_prime(s))-log(P(r))*2*N0/(N^2*(abs(H_f(s)))^2)];
            end
            if (s==r && norm(a(s)-a_prime(r))<1)               % If the distance is less than 1 we do not need to evaluate all 63 distances
                Certainly_NoError=1;
                Decoded=[Decoded a(s)];
                if figure_flag==1
                    figure(9)
                    plot(a(s),'b*')
                    hold on
                    plot(a_prime(r),'g*')
                    grid on
                    title('OFDM TX-RX constelation for ML receiver')
                    hold on
                end
                break
            end
        end
        if Certainly_NoError==0
            minimum_distance_index=find(distance_vector==min(distance_vector),1);
            Decoded=[Decoded SubCarrier_Mat(minimum_distance_index)];
            if SubCarrier_Mat(minimum_distance_index)==a(s)
                if figure_flag==1
                    figure(9)
                    plot(a(s),'b*')
                    hold on
                    plot(a_prime(s),'g*')
                    grid on
                    title('OFDM TX-RX constelation for ML receiver')
                    hold on
                end
            else
                if figure_flag==1
                    figure(9)
                    plot(a(s),'b*')
                    hold on
                    plot(a_prime(s),'r*')
                    grid on
                    title('OFDM TX-RX constelation for ML receiver')
                    hold on
                end
                if ML==1
                    NrSymbolErrorsML=NrSymbolErrorsML+1;
                end
                if MAP==1
                    NrSymbolErrorsMAP=NrSymbolErrorsMAP+1;
                end
            end
        end
    end

    % Decoding of method that involves message coding after each three     % interation
    if (Coding==1)
        A_Message(iterations).Received=a_prime;
        % % decoding of received message
        detected_message=[];
        for r=1:K
            detected_message=[detected_message find(abs(reshape(SubCarrier_Mat,64,1)-a_prime(r)).^2==min(abs(reshape(SubCarrier_Mat,64,1)-a_prime(r)).^2))];
        end
        A_Message(iterations).Detected_Vec=detected_message;
        A_Message(iterations).Detected_Message=reshape(detected_message,8,8);

    end
    if (Coding==1 &&  mod(iterations,3)==0)

        %we check at the receiver whether or not the
        %current symbol, for each subchannel, is the permuted linear combination of the previous two symbols
        %modulo 64 and addition with 1.


        ind_m1=iterations-2;  % information of first (main) message
        ind_m2=iterations-1;  % information of second (main) message
        ind_m3=iterations;    % information of third (redundant) message

        for r=1:K
            m1=A_Message(ind_m1).Detected_Message(r);
            m2=A_Message(ind_m2).Detected_Message(r);
            m3=A_Message(ind_m3).Detected_Message(r);

            % Chck the both linear combination and compare with m1 and m2
            if (m3==p(mod(m1+m2,64)+1) && m1==A_Message(ind_m1).Vec(r) && m2==A_Message(ind_m2).Vec(r))
                % there is no error :)
                stop_here_for_moment=1; % :)
            else
                c1_N=find(abs(A_Message(ind_m1).Received(r)-reshape(SubCarrier_Mat,1,64))<2*sqrt(2));  % Neighbors of (main) message 1
                c2_N=find(abs(A_Message(ind_m2).Received(r)-reshape(SubCarrier_Mat,1,64))<2*sqrt(2));  % Neighbors of (main) message 2
                c3_N=find(abs(A_Message(ind_m3).Received(r)-reshape(SubCarrier_Mat,1,64))<2*sqrt(2));  % Neighbors of (redundant) message 3

                Candidates=[];
                for c1=1:length(c1_N)
                    for c2=1:length(c2_N)
                        for c3=1:length(c3_N)
                            if c3_N(c3)==p(mod(c1_N(c1)+c2_N(c2),64)+1)
                                Candidates=[Candidates; c1_N(c1) c2_N(c2) c3_N(c3) abs(A_Message(ind_m1).Received(r)-SubCarrier_Mat(c1_N(c1))) abs(A_Message(ind_m2).Received(r)-SubCarrier_Mat(c2_N(c2))) abs(A_Message(ind_m3).Received(r)-SubCarrier_Mat(c3_N(c3))) abs(A_Message(ind_m1).Received(r)-SubCarrier_Mat(c1_N(c1)))+abs(A_Message(ind_m2).Received(r)-SubCarrier_Mat(c2_N(c2)))+abs(A_Message(ind_m3).Received(r)-SubCarrier_Mat(c3_N(c3)))];
                            end
                        end
                    end
                end
                Coding_Ind=find(Candidates(:,7)==min(Candidates(:,7)));   % finding the minimum distance satisfying coding conditions
                if (A_Message(ind_m1).Vec(r)~=Candidates(Coding_Ind,1))
                    NrSymbolErrorsCoding=NrSymbolErrorsCoding+1;
                end
                if (A_Message(ind_m2).Vec(r)~=Candidates(Coding_Ind,2))
                    NrSymbolErrorsCoding=NrSymbolErrorsCoding+1;
                end

            end
        end
    end

end

if Coding==1
    SymbolErrorsRate_Coding=NrSymbolErrorsCoding/(Nmes*K*2/3)
end
if ML==1
    SymbolErrorRate_ML=NrSymbolErrorsML/(Nmes*K)
end
if MAP==1
    SymbolErrorRate_MAP=NrSymbolErrorsMAP/(Nmes*K)
end





N0_Vec=[1000 1500 2000 2500 3000 3500 4000];
ML_Error=[0.00173 0.0111 0.0291 0.0541 0.0869 0.1202 0.151];
MAP_Error=[0.0023 0.0112 0.0303 0.0549 0.0827 0.116 0.147];
Coding_Error=[1.3e-4 1.95e-4 3.25e-4 0.0044 0.0108 0.0177 0.0304];
figure(10)
plot(N0_Vec,ML_Error,'-b*','LineWidth', 1.5)
hold on
plot(N0_Vec,MAP_Error,'-r*','LineWidth', 1.5)
grid on
plot(N0_Vec,Coding_Error,'-g*','LineWidth', 1.5)
grid on
legend('ML','MAP','Coding')



