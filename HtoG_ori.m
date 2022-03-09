keep_going = 0;
 
    H = matrix;
 
    column_weights = sum(H,1);
 
    H2 = H;
   [m,n] = size(H2);
   
   nminusk = m;
   k = n-m;
   
    for row_index_1 = 1:nminusk
        % Find the first row which can be swapped to remove a zero on the diagonal
        row_index_2 = row_index_1;
        while row_index_2 <= nminusk && H2(row_index_2, row_index_1) == 0
            row_index_2 = row_index_2 + 1;
        end
        
        if row_index_2 <= nminusk
        
            % Swap the rows if necessary
            if row_index_1 ~= row_index_2
                H2([row_index_1, row_index_2],:) = H2([row_index_2, row_index_1],:);
            end
 
        else
            column_index = row_index_1;
            while column_index <= n && ~(H2(row_index_1, column_index) == 1 && column_weights(row_index_1) == column_weights(column_index))
                column_index = column_index + 1;
            end
            
            if column_index <= n
                H2(:,[row_index_1, column_index]) = H2(:,[column_index, row_index_1]);
                H(:,[row_index_1, column_index]) = H(:,[column_index, row_index_1]);
            else
                keep_going = 1;
                row_index_1 = nminusk + 1;
            end
            
        end
        
        if keep_going == 0
            % Convert the rows below into upper triangular form
            for row_index_2 = 1:row_index_1-1
                if H2(row_index_2,row_index_1)
                    H2(row_index_2,:) = mod(H2(row_index_2,:) + H2(row_index_1,:),2);
                end
            end
            for row_index_2 = row_index_1+1:nminusk
                if H2(row_index_2,row_index_1)
                    H2(row_index_2,:) = mod(H2(row_index_2,:) + H2(row_index_1,:),2);
                end
            end
        end
            
    end    
    P = H2(:,nminusk+1:end);
    %H3 = [P,eye(m)]
    G = [P',eye(k)];
    
    
      [r,k] = size(G);
      
    for attempt = 1:10000  
      
    X = randi([0 1],1,r); 
    C = mod(X*G,2);
    if ~isempty(find(mod(C*H.',2)~=0))
        error('error');
    end
     end
