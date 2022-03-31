function P_proj = Simplex_proj(P,S)

% A function which returns the projection onto a simplex
% Projection Onto A Simplex, Yunmei Chen and Xiaojing Ye, 2011


  
    t_vec = zeros(S-1);

    % Step 2
    for i = 1:S
        P_sort(i,:) = sort(P(i,:));    
    end
    P_sort;
    % Step 3
    for i = 1:S
        for j = 1:S-1
            t_vec(i,j) = (sum(P_sort(i,j+1:end))-1)/(S-j);
            (sum(P_sort(i,j+1:end))-1)/(S-j)
        end
    end
    
    % t_vec
    
    % P_sort
    
    for i = 1:S
        flag = zeros(S-1);
        for j_temp = 1:S-1
            if t_vec(i,S-j_temp)>=P_sort(i,S-j_temp)          
                flag(S - j_temp) = 1;
            else               
            end
        end
        if sum(flag) == 0
            hat_t(i) = sum(P_sort(i,:)-1)/S;
        else
            pos = find(flag);
            hat_t(i) = t_vec(i,pos(end));
        end
        P_proj(i,:) = max(P(i,:) - hat_t(i) * ones(1,S),0);
    end
    hat_t;
%     for i = 1:S
%         for j_temp = 1-S
%         max (t_vec(i,:)>=P_sort(i,S-j_temp))
%                 hat_t(i) = t_vec(i,S - j_temp);
%                 j_temp = S; % a local break
%             else               
%             end
%         end
%         if j_temp == S
%             hat_t(i) = sum(P_sort(i,:)-1)/S;
%         else
%         end
%         P_proj(i,:) = max(P(i,:) - hat_t(i) * ones(1,S),0);
%     end
    
end