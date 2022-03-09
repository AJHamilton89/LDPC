function [codewordout,Codeword_correct, syndrome] = LDPC_Decoder_Standalone(Channel,matrix,num_iterations)

H=matrix;
num_v = size(H,2);
num_c = size(H,1);

Codeword_correct=0;


   
%         tic
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    
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
            mutual_information_c(n) = 1-sum(entropies_c(find (matrix > 0)))/nnz(matrix);
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
    
    
%     P0 = exp(Channel)./(1+exp(Channel));
%     P1 = 1-P0;
%     entropies = -P0.*log2(P0)-P1.*log2(P1);
%     MI_Channel = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
      
    %% display LLRS
%       disp(sprintf('LLRs after %d iterations: %s',n,mat2str(sums,2)));
%     processedLLRstore(Noise_lvls,l)=sum(abs(sums));
    
    %% display EXIT
%     if n>=1
%     if (mod(xhat*matrix.',2)) == 0
%     MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v 1];
%     MI_EXIT_C = [0 mutual_information_c 1];
%     else
%         MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v];
%     MI_EXIT_C = [0 mutual_information_c];
%     end
%     figure (1)
%     cla
%     hold off
%     axis equal
%     grid minor
%     xlim([0 1])
%     ylim([0 1])
%     title('Full EXIT Function')
%     xlabel('I_A')
%     ylabel('I_E')
%     hold on
%     plot (MI_EXIT_C,MI_EXIT_V,'d--'),stairs(MI_EXIT_C,MI_EXIT_V,'-')
%     
%     %find equivaelnt CND EXIT function
%     [xCND,yCND]=stairs(MI_EXIT_C,MI_EXIT_V);
%     MI_CNDx=xCND(2:2:end);
%      MI_CNDy=yCND(2:2:end);
%      plot(MI_CNDx,MI_CNDy,'+-.')
%      
%      drawnow
%     
%     end
    %% Calculate Syndrome
     
    s = (mod(xhat*matrix.',2));
    syndrome=s;
    codewordout=xhat;
    
    if (sum(s) == 0)
%     disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
    Codeword_correct=1;
    end
%     if (sum(s) == 0)
    
        
        
%         %[~,indx] = ismember(xhat,codewords,'rows');
%     %if (indx > 0)
%         disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
%         correctcount = correctcount + 1;
%         if sum(xhat-codeword)~=0
%             falsecodeword=falsecodeword+1
%         end
%          if sum(zhat-xhat)==0
%             LDPCuseless=LDPCuseless+1
%         end
%     else
%         disp(['WRONG CODEWORD: ' mat2str(xhat) sprintf('\n')]);
    


    end