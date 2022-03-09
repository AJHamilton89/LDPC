% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

% QAM - 4 gives QPSK, 64 gives 64QAM
% R - see puncture.m for supported coding rates
% K - see puncture.m for supported info block lengths
% channel - 'AWGN' gives an AWGN channel, 'Rayleigh' gives an uncorrelated narrowband Rayleigh channel
% EsN0_dB_start - first EsN0 value to use
% EsN0_dB_delta - increment to use for subsequent EsN0 values
% new_approx_maxstar -  0 gives Log-MAP, 1 gives Max-Log-MAP, 2 gives scaled-Max-Log-MAP, 3 gives constant-Max-Log-MAP, 4 gives scaled-constant-Max-Log-MAP
% max_iterations - number of turbo decoder iterations to perform
% block_errors_limit - simulation continues until this many errors have been observed
function main_bler_turbo(QAM,R,K,channel,EsN0_dB_start,EsN0_dB_delta,new_approx_maxstar,max_iterations,block_errors_limit)

if nargin < 1
    QAM = 4;
end
if nargin < 2
    R = 0.2;
end
if nargin < 3
    K = 992;
end
if nargin < 4
    channel = 'AWGN';
end
if nargin < 5
    EsN0_dB_start = -5; % Choose the SNR
end
if nargin < 6
    EsN0_dB_delta = 0.5;
end
if nargin < 7
    new_approx_maxstar = 2;
end
if nargin < 8
    max_iterations = 10;
end
if nargin < 9
    block_errors_limit = 100;
end

global approx_maxstar
approx_maxstar=new_approx_maxstar;
iterations = 0:max_iterations;



filename = sprintf('results_turbo_%d_%f_%d_%s_%f_%f_%d_%d_%d.txt',QAM,R,K,channel,EsN0_dB_start,EsN0_dB_delta,new_approx_maxstar,max_iterations,block_errors_limit);
fid = fopen(filename,'w');
fprintf(fid,'QAM\tRate\tInfo Blocklength\tEs/N0\tBLER\tChannel\n');


% Create a figure to plot the results.
figure
axes1 = axes('YScale','log');
    title(['Enhanced Turbo ',num2str(QAM),'QAM R = ',num2str(R),' K = ',num2str(K), ' ', channel]);
ylabel('BLER');
xlabel('E_s/N_0 [dB]');
hold on
plots = zeros(size(iterations));
for iteration_index = 1:length(iterations)
    plots(iteration_index) = plot(nan,'Parent',axes1);
end

EsN0_dB = EsN0_dB_start;
block_errors = zeros(1,length(iterations));
blocks = 0;

while 1
    
matrix = genmatrix('rb_alpha.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;
%
int_LLRs=zeros(size(H));

correctcount = 0;
falsecodeword=0;
LDPCuseless=0;
num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
syndrome_output=10;
num_runs=1;

% set scaling
minSNR=0;
SNRstepsize=1;
num_points=5;
maxSNR=minSNR+SNRstepsize*(num_points-1);

[r,k] = size(G);
for Noise_lvls=1:num_points
    for l=1:num_runs
%% Generate Codeword using G Matrix (Note:13/1/2017 unsure if G is correct)
 X = randi([0 1],1,r); 
 codeword = mod(X*G,2);
 %codeword = DFENCRB;
 
    %% Generate random noise variance in accordance with Channel model and modulation
    
SNR = minSNR+((Noise_lvls-1)*SNRstepsize)
     N0 = 1/(10^(SNR/10));

     a_tx = 2*(codeword-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+i*randn(size(a_tx)));
     Channel = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;
    
    disp(['codeword: ' mat2str(codeword)]);
    disp(['received LLRs: ' mat2str(Channel,2)]);
    
    %% decision on unprocessed codeword based on LLRS 
    
    rawLLRstore(Noise_lvls,l)=sum(abs(Channel));
    zhat = zeros(1,num_v);
    zhat(Channel>0) = 1;

    %% Load LLRs into Parity Check
    H(H==0) = NaN;
    for c = 1:num_c
        for v = 1:num_v
            if matrix(c,v)==1
                H(c,v) = Channel(v);
            end
        end
    end
    
  
    %% start iterative decoding
    H(isnan(H))=0;
    clear MI_v
    clear mutual_information_c
    
    for n = 1:num_iterations;
%         tic
            HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    
    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums>0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
  end
        %% check node calculations
        for c = 1:num_c
            int_LLRs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v
          
                if matrix(c,v)==1
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
        
        P0_c = exp(H)./(1+exp(H));
            P1_c = 1-P0_c;
            entropies_c = -P0_c.*log2(P0_c)-P1_c.*log2(P1_c);
            mutual_information_c(n) = 1-sum(entropies_c(~isnan(H)))/numel(entropies_c);
        %% variable node calculations
       
        for v = 1:num_v
            int_LLRs = Channel(v);
            for c = 1:num_c
                if matrix(c,v)==1
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = variable_node(int_LLRs);
            
           % pri_LLRs
            pri_LLRs(1) = []; %remove channel value before re-writing to H
           
            for c = 1:num_c
                
                if matrix(c,v)==1
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
        %disp(['end of iteration ',num2str(n)]);
%         toc
        H(isnan(H)) = 0;
        
        P0_v = exp(nonzeros(H))./(1+exp(nonzeros(H)));
            P1_v = 1-P0_v;
            entropies_v = -P0_v.*log2(P0_v)-P1_v.*log2(P1_v);
            MI_v(n) = 1-sum(sum(entropies_v(~isnan(entropies_v)),1))/numel(entropies_v);
  HplusC = [Channel;H];
    sums = sum(HplusC);
    end
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
       xhat = zeros(1,num_v);
    xhat(sums>0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
   
    %% what is the Mutual information at the end of the iteration
    P0 = exp(sums)./(1+exp(sums));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    MI_c(Noise_lvls,n) = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
    
    P0 = exp(Channel)./(1+exp(Channel));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    MI_Channel(Noise_lvls,l) = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
      
    %% display LLRS
      disp(sprintf('LLRs after %d iterations: %s',n,mat2str(sums,2)));
    processedLLRstore(Noise_lvls,l)=sum(abs(sums));
    
    %% display EXIT
    if n<=1
    if (mod(xhat*matrix.',2)) == 0
    MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v 1]
    MI_EXIT_C = [0 mutual_information_c 1]
    else
        MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v]
    MI_EXIT_C = [0 mutual_information_c]
    end
    figure
    axis equal
    grid minor
    xlim([0 1])
    ylim([0 1])
    title('Full EXIT Function')
    xlabel('I_A')
    ylabel('I_E')
    hold on
    plot (MI_EXIT_C,MI_EXIT_V,'d--'),stairs(MI_EXIT_C,MI_EXIT_V,'-')
    
    %find equivaelnt CND EXIT function
    [xCND,yCND]=stairs(MI_EXIT_C,MI_EXIT_V);
    MI_CNDx=xCND(2:2:end);
     MI_CNDy=yCND(2:2:end)
     plot(MI_CNDx,MI_CNDy,'+-.')
    end
    %% Calculate Syndrome
     
    s = (mod(xhat*matrix.',2));
    if (sum(s) == 0)
    %[~,indx] = ismember(xhat,codewords,'rows');
    %if (indx > 0)
        disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
        correctcount = correctcount + 1;
        if sum(xhat-codeword)~=0
            falsecodeword=falsecodeword+1
        end
         if sum(zhat-xhat)==0
            LDPCuseless=LDPCuseless+1
            
                        block_errors(end) = block_errors(end) + 1;
                   
        end
    else
        disp(['WRONG CODEWORD: ' mat2str(xhat) sprintf('\n')]);
    


    end
    correctcount
    
    %% Gather Data at end of iterations
    rawErrors(Noise_lvls,l)=sum(abs(zhat-codeword));
    Errors(Noise_lvls,l)=sum(abs(xhat-codeword));
    Syndrome(Noise_lvls,l)=sum(s);
    
    end
end
for iteration_index = 1:length(iterations)
            set(plots(iteration_index),'XData',EsN0_dB);
            set(plots(iteration_index),'YData',block_errors(:,iteration_index)./blocks);
        end
        drawnow
figure
end

% scale = minSNR:SNRstepsize:maxSNR;
% rerr=sum(rawErrors,2);
% rBER=rerr/(r*l);
% semilogy(scale,rBER)
% title('Raw Bit Error Rate in AWGN Channel')
% xlabel('SNR')
% ylabel('Raw Bit Error Rate')



% figure
% 
% err=sum(Errors,2);
% BER=err/(r*l);
% semilogy(scale,BER)
% title('Bit Error Rate for AWGN Channel after LDPC')
% xlabel('SNR')
% ylabel('Bit Error Rate')
% 
% figure
% 
% syn=sum(Syndrome,2);
% SER=syn/(r*l);
% semilogy(scale,SER)
% title('Syndrome Error Rate for AWGN Channel')
% xlabel('SNR')
% ylabel('Syndrome Error Rate')
% 
% 
% figure
% 
% plot(Errors,Syndrome,'.')
% title('Correlation between Bit Errors and Syndrome Errors in AWGN channel w/BPSK')
% xlabel('Bit Error Rate')
% ylabel('Syndrome Error Rate')
% 
% figure
% 
% plot(BER,SER,'.')
% title('Correlation between Bit Errors and Syndrome Error Rates in AWGN channel w/BPSK (averaged)')
% xlabel('Bit Error Rate')
% ylabel('Syndrome Error Rate')
% 
% figure
% 
% plot(scale,MI_c)
% title('Mutual Information of LLRs post LDPC')
% xlabel('SNR')
% ylabel('Final MI')
% 
% 
% figure
% 
% plot(scale,MI_Channel)
% title('Mutual Information of LLRs at output of Modulator')
% xlabel('SNR')
% ylabel('Intial MI')
% 
% 
% figure
% 
% plot(MI_Channel,MI_c,'.')
% title('Relationship between MIs of Channel LLRs and post-LDPC LLRs')
% xlabel('Channel MI')
% ylabel('Output MI')
