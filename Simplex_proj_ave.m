function t = Simplex_proj_ave(p,i)

% Projection Onto A Simplex, Yunmei Chen and Xiaojing Ye, 2011
% Step 3 

    Size_p = size(p);
    S = Size_p(1);

    t(i) = (sum(P_sort(i+1:end,:))-1)/(S-i);

end
