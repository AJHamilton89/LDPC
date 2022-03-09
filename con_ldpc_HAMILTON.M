
matrix = genmatrix('rb_bravo.txt');
[G,H3] = HtoG(matrix);
%H=matrix;%
H = H3;
%
int_LLRs=zeros(size(H));
matrix=H3;
correctcount = 0;
num_v = size(H,2);
num_c = size(H,1);
num_iterations = 10;
syndrome_output=10;

[r,k] = size(G);
for SNR=6:8
    for l=1:1
%% Generate Codeword using G Matrix (Note:13/1/2017 unsure if G is correct)
 X = randi([0 1],1,r); 
 codeword = mod(X*G,2);
 
    %% Generate random noise variance in accordance with Channel model and modulation
    

     N0 = 1/(10^(SNR/10));

     a_tx = 2*(codeword-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+i*randn(size(a_tx)));
     Channel = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;
    
    disp(['codeword: ' mat2str(codeword)]);
    disp(['received LLRs: ' mat2str(Channel,2)]);
    
    %% decision on unprocessed codeword based on LLRS 
    
    rawLLRstore(SNR,l)=sum(abs(Channel));
    zhat = zeros(1,num_v);
    zhat(Channel>0) = 1;

    %% Load LLRs into Parity Check
    H(H==0) = NaN;
    for c = 1:num_c
        for v = 1:num_v
            if (~isnan(H(c,v)))
                H(c,v) = Channel(v);
            end
        end
    end
    
    %% start iterative decoding
    H(isnan(H))=0;
    
    for n = 1:num_iterations;
        tic
            HplusC = [Channel;H];
    sums = sum(HplusC);
    
    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums>0) = 1;
    s = (mod(xhat*H3.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
  end
        %% check node calculations
        for c = 1:num_c
            int_LLRs = [];
            for v = 1:num_v
                if (~isnan(H(c,v)))
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
           % pri_LLRs
            for v = 1:num_v
                if (~isnan(H(c,v)))
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        end
        
        %% variable node calculations
       
        for v = 1:num_v
            int_LLRs = Channel(v);
            for c = 1:num_c
                if (~isnan(H(c,v)))
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = variable_node(int_LLRs);
           % pri_LLRs
            pri_LLRs(1) = []; %remove channel value before re-writing to H
            for c = 1:num_c
                if (~isnan(H(c,v)))
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        end
        %disp(['end of iteration ',num2str(n)]);
        toc
        H(isnan(H)) = 0;

    end
    
    
      
    %% display LLRS
      disp(sprintf('LLRs after %d iterations: %s',n,mat2str(sums,2)));
    processedLLRstore(SNR,l)=sum(abs(sums));
    
    
    
    %% Calculate Syndrome
    
    s = (mod(xhat*H3.',2));
    if (sum(s) == 0)
    %[~,indx] = ismember(xhat,codewords,'rows');
    %if (indx > 0)
        disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
        correctcount = correctcount + 1;
    else
        disp(['WRONG CODEWORD: ' mat2str(xhat) sprintf('\n')]);
    


    end
    correctcount
    
    %% Gather Data at end of iterations
    rawErrors(SNR,l)=sum(abs(zhat-codeword));
    Errors(SNR,l)=sum(abs(xhat-codeword));
    Syndrome(SNR,l)=sum(s);
    end
end
figure

rerr=sum(rawErrors,2);
rBER=rerr/(r*l);
semilogy(rBER)
title('Raw Bit Error Rate in AWGN Channel')
xlabel('SNR')
ylabel('Raw Bit Error Rate')

figure

err=sum(Errors,2);
BER=err/(r*l);
semilogy(BER)
title('Bit Error Rate for AWGN Channel after 20 iterations')
xlabel('SNR')
ylabel('Bit Error Rate')

figure

syn=sum(Syndrome,2);
SER=syn/(r*l);
semilogy(SER)
title('Syndrome Error Rate for AWGN Channel after 20 iterations')
xlabel('SNR')
ylabel('Syndrome Error Rate')

figure

plot(BER,SER'.')
title('Correlation between Bit Errors and Syndrome Errors in AWGN channel w/BPSK')
xlabel('Bit Error Rate')
ylabel('Syndrome Error Rate')
