function ber

    %clc;
    %close all;
    %clear all;
%% LDPC Initilastion functions


matrix = genmatrix('rb_alpha.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;
%Moll

for number_iterations = 1:10:50:100:200

[r,k] = size(G);


%% rest of code



    
    %omp = 1;
    %carrier = 0;
    
    SNR_range = 0:1:20;
    snr_linear=10.^(SNR_range/10);
    n_frames = 1000;
    
    n_errors_perfectcsi = 0;
    total_bits_perfectcsi = 0;
    ber1 = [];
    
    n_bl_er_perf = 0;
    bler1 = [];
    
    T=2*4/3*0.384; % time of the signal
    fft_size = 64;
    cyclicPrefix = fft_size/8;
    %pilotDensity = 0.2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set the system parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_s = 64000;    % sampling frequency
    f_c = 25000;    % carrier frequency
    f_b = 8000;     % bit rate
    f_b = k; % EDITED - default bit rate to LDPC packet size
    
    % set the number of bits for the simulation
    bits=(lcm(fft_size,k/2));
%     bits = f_b * T;
    b_hat = [];
    
    pskModulator = comm.PSKModulator(4,-pi/4,'BitInput',true);
    %pskDemodulator = comm.PSKDemodulator(4,-pi/4,'BitOutput',true);
    pskDemodulator = comm.PSKDemodulator('ModulationOrder',4,...
        'PhaseOffset',-pi/4,'BitOutput',true, 'DecisionMethod', 'Log-likelihood ratio');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BER/SNR simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %snr_i = 1;
    for snr_i = 1:length(SNR_range)
        SNR = SNR_range(snr_i)
        %snr_i = snr_i + 1;
        N0 = 1/snr_linear(snr_i);
        
        n_errors_perfectcsi = 0;
        total_bits_perfectcsi = 0;
        n_bl_er_perf = 0;
        
        for frames = 1:n_frames
            
            % generate bits
            b=[]; 
            for n=1:(lcm(fft_size,k/2)/(k/2))
                 
                 X = randi([0 1],1,r); 
            codeword = mod(X*G,2); %LDPC codeword
            b= [b codeword];
            end
            
            %b = round(rand(1,bits)); %rounded bits
            
            % modulator specs
            bits_per_symbol = 2;                % QPSK
            symbols = bits/bits_per_symbol;     % number of symbols
            f_symb = f_b/bits_per_symbol;       % symbol rate
            samples_per_symbol = f_s/f_symb;    % number of samples per symbols, which depends on the sampling frequency
            
            samples=f_s*T; % number of samples
            t=(0:(samples-1))/f_s; % samples times
            
            % modulate the bits
%             b_mat = reshape(b, bits_per_symbol, symbols);
%             b_1 = b_mat(1,:);
%             b_2 = b_mat(2,:);
%             a=sqrt(1/2);
%             sn = a*(-2*(b_1-0.5) + 1j*-2*(b_2-0.5));
            
            sn = pskModulator(b.');
            
            % form the OFDM symbol
            CP = cyclicPrefix;
            
            % number of OFDM symbols
            Nsyms = length(sn)/fft_size;
            % reshape to form the OFDM symbols
            sn1 = reshape(sn, fft_size, Nsyms).';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ber_perfectCSI = 0;
            b_1hat = [];
            b_2hat = [];
            b_hat = [];
            channel_llrs = [];
            
            for ab = 1:size(sn1,1)
                % Taking IFFT
                x_ifft(ab,:) = sqrt(fft_size) * ifft((sn1(ab,:)),fft_size);
                
                % add cyclic prefix
                if (CP <= fft_size)
                    x_cp = [x_ifft(ab,[end-CP+1:end]) x_ifft(ab,:)];
                    x_cp1 = x_cp;
                else
                    nn = ceil(CP / fft_size);
                    x_cp = [x_ifft(ab,:)];
                    for io = 1:nn
                        x_cp = [x_cp x_ifft(ab,:)];
                    end
                    x_cp1 = x_cp;
                end
                
                  % define UWA channel
                % set nzt non-zero taps out of notap taps
                nzt=1;
                notap=1;
                
                pp1 = ceil(notap*rand(1,nzt)); 
                hh1 = (1/sqrt(2))* (randn(1,nzt) + 1j*randn(1,nzt));
                hh1 = hh1 / sqrt(sum(abs(hh1)));
                h = zeros(1,notap);
                h(pp1) = hh1;
                channelTD = h;
                                
                % transmit through the channel
                yhat1 = conv(x_cp1,channelTD);
                
%                 figure(1)
%                 plot(abs(channelTD))
%                 
                % add noise
                n =  sqrt(N0/2)*(randn(1,length(yhat1))+1j*randn(1,length(yhat1)));
                yhat = yhat1 + n;
                
                % Remove cyclic prefix
                yc = yhat(:,[CP+1:size(x_ifft,2)+CP]);
                
                % transform to frequency domain
                y_fft(ab,:) = (1/sqrt(fft_size))*(fft(yc,fft_size));
                yc = y_fft(ab,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Get the perfect channel knowledge after the filtering %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cfd(ab,:) = (1/sqrt(fft_size))*(fft(channelTD,fft_size));
                ctd(ab,:) = sqrt(fft_size) * ifft(cfd(ab,:), fft_size);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % equalise using perfect CSI %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                yEq = yc ./ cfd(ab,:);
                % demodulate the signal
%                 b_1hat = [b_1hat (real(yEq)<0)];
%                 b_2hat = [b_2hat (imag(yEq)<0)];

                %b_hat = [b_hat; pskDemodulator(yEq.')];
                channel=pskDemodulator(yEq.') ;
                channel_llrs = [channel_llrs; channel];
                
                rxData = pskDemodulator(yEq.') < 0;
                b_hat = [b_hat; rxData];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % evaluate the BER using the perfect CSI %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            b_mathat = [b_1hat;b_2hat];
%            b_hat=reshape(b_mathat,1,2*length(b_1hat));
%            find(b_hat.' ~= b(1:length(b_hat)));
            n_errors_perfectcsi = n_errors_perfectcsi + sum(b_hat.' ~= b(1:length(b_hat)));
            total_bits_perfectcsi = total_bits_perfectcsi + length(b_hat);
            
            if (sum(b_hat.' ~= b(1:length(b_hat))) > 0)
                
                %% decoding of LDCP goes here need to get access to LLRs and feed into LDPC/Polar decoder
                channel_llrs=reshape(channel_llrs,[],size(channel_llrs,1)/k);
                
               
                for ldpcframe=1:size(channel_llrs,1)/k
                [outputbits,Codeword_correct,syndrome]=LDPC_Decoder_Standalone(channel_llrs(:,ldpcframe).',matrix,number_iterations);
                if Codeword_correct==0
                    n_bl_er_perf = n_bl_er_perf + 1;
                end
                end
                
%                 n_bl_er_perf = n_bl_er_perf + 1;
            else
                aaaa=1;
            end
            
        end
        ber1 = [ber1 n_errors_perfectcsi/total_bits_perfectcsi];
        
        % there are size(channel_llrs,1)/k (32) LDPC frames within one OFDM
        % frame
        bler1 = [bler1 n_bl_er_perf/(n_frames*size(channel_llrs,1)/k)];
    end
    
%    save('ber.mat','ber1','ber2', 'ber3', 'bler1','bler2', 'bler3')
    
figure (5)
hold on
    ber_perfCsi = ber1;
    semilogy(SNR_range,ber1)
    legend('OFDM system')
    xlabel('SNR [dB]')
    ylabel('BER')
    hold off
 
    
    figure(6);
    hold on
    semilogy(SNR_range, bler1)
    legend([num2str(number_iterations)])
    xlabel('SNR [dB]')
    ylabel('BLER')
    hold off
    
end
end
