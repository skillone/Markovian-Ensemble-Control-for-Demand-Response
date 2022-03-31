function P_proj = Simplex_proj_vec(P,S)

% A function which returns the projection onto a simplex
% Projection Onto A Simplex, Yunmei Chen and Xiaojing Ye, 2011


t_vec = zeros(1,S-1);

% Step 2
P_sort = sort(P(:));

P_sort;
% Step 3
for j = 1:S-1
    t_vec(j) = (sum(P_sort(j+1:end))-1)/(S-j);
end
t_vec;
flag = zeros(S-1);
for j_temp = 1:S-1
    if t_vec(S-j_temp)>=P_sort(S-j_temp)
        flag(S - j_temp) = 1;
    else
    end
end
if sum(flag) == 0
    hat_t = (sum(P_sort(:))-1)/S;
else
    pos = find(flag);
    hat_t = t_vec(pos(end));
end
P_proj = max(P - hat_t * ones(1,S),0);

P_proj;
end
%     for i = 1:S
%         for j_temp = 1-S
%         max (t_vec(:)>=P_sort(S-j_temp))
%                 hat_t(i) = t_vec(S - j_temp);
%                 j_temp = S; % a local break
%             else
%             end
%         end
%         if j_temp == S
%             hat_t(i) = sum(P_sort(:)-1)/S;
%         else
%         end
%         P_proj(:) = max(P(:) - hat_t(i) * ones(1,S),0);
%     end

