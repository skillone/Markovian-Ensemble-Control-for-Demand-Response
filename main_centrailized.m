
%% Solve for value functions of each unit

% Ensemble state vector
x_default = zeros(S_num, L, N);
x_controlled = zeros(S_num, L, N);

% Value function
global V
V = zeros(S_num,L-1, N);
% Some intermediate value
z_temp = zeros(S_num, L-1, N);
for n = 1:N
    
    % Calculate Z value (Eq.3 from paper)
    z = zeros(S_num, L);
    z(:, L) = exp(ele_cost(:, L) ./ gamma(n));
    for l = 1: L-1
        for i = 1:S_num
            temp1 = 0;
            for j = 1:S_num
                temp1 = temp1 + z(j, L-l+1) * P_default_N(i, j, n); % Summation in Eq.3
            end
            z(i, L-l) = temp1 * exp(-ele_cost(i, L-l) / gamma(n)); % Final Eq.3
        end
    end
    
    % Compute for value function
    for l = 1:L-1
        for i = 1:S_num
            for j = 1:S_num
                z_temp(i,l,n) = P_default_N(i, j, n)*z(j, l)+z_temp(i,l,n);
            end
        end
    end
    
    for l = 1:L-1
        for i = 1:S_num
            V(i,l,n) = ele_cost(i,l) -gamma(n) * log(z_temp(i,l,n));
        end
    end
    V(:,L,n) = ele_cost(:,L);


    
end






%% Centralized Scheme

gamma_bar = 0;
    for n = 1:N
        gamma_bar = gamma(n) + gamma_bar;
    end
gamma_bar = gamma_bar./N;


mu_bar = zeros(S_num,S_num,N);
for i = 1:S_num
    for j = 1:S_num
        for n = 1:N
            mu_bar(i,j) = gamma(n)*log(P_default_N(i,j,n)) + mu_bar(i,j);
        end
    end
end
mu_bar = mu_bar./N;


% Calculate Z value (Eq.3 from paper)
v_bar = zeros(S_num, L);
v_bar(:, L) = ele_cost(:,L);

V_ave = zeros(S_num,L);
V_ave(:, L) = ele_cost(:,L); 
for l = 1:L-1
    for i = 1:S_num
        V_ave(i,l) = sum(V(i,l,:))/N;
    end
end

for l = 1: L-1
    for i = 1:S_num
        temp1 = 0;
        for j = 1:S_num
            temp1 = temp1 + exp((mu_bar(i,j)-V_ave(j, L-l+1))/gamma_bar); % Summation in Eq.3
        end
        v_bar(i,L-l) = ele_cost(i,L-l) - gamma_bar * log(temp1); % Final Eq.3
    end
end

V_ave = zeros(S_num,L);
V_ave(:, L) = ele_cost(:,L); 
for l = 1:L-1
    for i = 1:S_num
        V_ave(i,l) = sum(V(i,l,:))/N;
    end
end


%Calculate P (Eq.2 from paper)
P_control_bar = zeros(S_num, S_num, L - 1);
for l = 1:L-1
    for i = 1:(S_num)
        temp2 = 0;
        for j = 1:(S_num)
            temp2 = temp2 +  exp((mu_bar(i,j)-v_bar(j, l+1))/gamma_bar); %denominator in Eq.2
        end
        for j = 1:(S_num)
            P_control_bar(i, j, l) = exp((mu_bar(i,j)-v_bar(j, l+1))/gamma_bar) / temp2; %Final Eq.2
        end
    end
end


% calculating x (default and controlled) (probability of each state at time t)
x_default_bar = zeros(S_num, L);
x_controlled_bar = zeros(S_num, L);
x_default_bar(:, 1)= P_default(S_num, :); % initial state
x_controlled_bar(:, 1)= P_default(S_num, :); % initial state

% Default ensemble vector
for l = 1:(L-1)
    for i = 1:S_num
        temp3 = 0;
        for j = 1:S_num
            temp3 = temp3 + P_default_N(j, i) * x_default_bar(j, l);
        end
        x_default_bar(i, l+1) = temp3;
    end
end

% Controlled ensemble vector
for l = 1:(L-1)
    for i = 1:S_num
        temp3 = 0;
        for j = 1:S_num
            temp3 = temp3 + P_control_bar(j, i, l) * x_controlled_bar(j, l);
        end
        x_controlled_bar(i, l+1) = temp3;
    end
end

% Electricity consumption (default and controlled)
for l = 1:L
    elec_use_default_bar(l) = elec_use * x_default_bar(:,l);
    elec_use_controlled_bar(l) = elec_use * x_controlled_bar(:,l);
end

save('elec_use_controlled_bar')


