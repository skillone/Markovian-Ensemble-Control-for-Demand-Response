%% Solve for value functions of each unit

% Ensemble state vector
x_default = zeros(S_num, L, N);
x_controlled = zeros(S_num, L, N);


%% Generate communication matrix
G = ones(N,N)*(1/N);

% Generate random default matrix
beta_commm = 0.8; % parameter 1-beta is the level of how much the units are subject to randomness
for n = 1:N
    for m = 1:N
        r_comm(m,:) = rand(1,N); % generate random row vector
        r_comm_temp(m) = sum(r_comm(m,:));  % Normalize so the sum is 1.
    end
    
    for m = 1:N
        for j = 1:S_num
            %             G_mat(m,n) = G(m,n) * beta_commm + (1-beta_commm) * r_comm(m,n)/r_comm_temp(m);
            G_mat(m,n) =  1/N;
        end
    end
end

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
            V(i,l,n) = ele_cost(i,l) - gamma(n) * log(z_temp(i,l,n));
        end
    end
    V(:,L,n) = ele_cost(:,L);
    
    % Calculate P (Eq.2 from paper)
    P_control = zeros(S_num, S_num, L - 1);
    for l = 1:L-1
        for i = 1:(S_num)
            temp2 = 0;
            for t = 1:(S_num)
                temp2 = temp2 + z(t,L-l+1) * P_default_N(i,t,n); %denominator in Eq.2
            end
            for j = 1:(S_num)
                P_control(i,j,L-l) = z(j,L-l+1) * P_default_N(i,j,n) / temp2; %Final Eq.2
            end
        end
    end
    
    
    % calculating x (default and controlled) (probability of each state at time t)
    x_default(:, 1, n)= P_default_N(S_num, :, n); % initial state
    x_controlled(:, 1, n)= P_default_N(S_num, :, n); % initial state
    
    % Default ensemble vector
    for l = 1:(L-1)
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + P_default_N(j, i, n) * x_default(j, l, n);
            end
            x_default(i, l+1, n) = temp3;
        end
    end
    
    % Controlled ensemble vector
    for l = 1:(L-1)
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + P_control(j, i, l) * x_controlled(j, l, n);
            end
            x_controlled(i,l+1,n) = temp3;
        end
    end
    
    % Electricity consumption (default and controlled)
    for l = 1:L
        elec_use_default(l,n) = elec_use * x_default(:,l,n);
        elec_use_controlled(l,n) = elec_use * x_controlled(:,l,n);
    end
    
    
    % store the optimal control at a given stage
    for l = 1 :L-1
        P_control_vec = reshape(P_control(:, :, l)',[S_num*S_num,1]);
        Vec_State(:,l,n) = P_control_vec;
    end
    
end

P_default_ave = zeros(S_num,S_num);
for n = 1:N
    P_default_ave(:, :) = P_default_ave(:, :) + P_default_N(:, :, n);
end
P_default_ave = P_default_ave./N;
for i = 1:S_num
    for j = 1:S_num
        if P_default_ave(i, j) == 0
            P_default_ave(i, j) = epsilon; % a small probability
        end
    end
end

x_default_ave(:, 1)= P_default_ave(S_num, :); % initial state
% Controlled ensemble vector
for l = 1:(L-1)
    for i = 1:S_num
        temp3 = 0;
        for j = 1:S_num
            temp3 = temp3 + P_default_ave(j, i) * x_default_ave(j,l);
        end
        x_default_ave(i,l+1) = temp3;
    end
end
for l = 1:L
    elec_use_default_ave(l) = elec_use * x_default_ave(:,l);
end

% 
% % Invidual controlled
% figure(2)
% time = 1:L;
% for n = 1 : N
%     plot(time, elec_use_default(:,n), '-k')
%     hold on
%     plot(time, elec_use_controlled(:,n), '-.r')
% end
% %legend('PEC Consumption','Observed Consumption','Prior-Estimated Consumption','Posterior-Estimated Consumption')
% plot(time, elec_use_default_ave, '-b', 'LineWidth',2)
% set(gca,'FontSize',15)
% grid on
% xlabel('Time (min)')
% ylabel('Power (MW)')


%% Decentralized Setting

% A test for local consensus

T = 10000; % time index


t = 1;


err_Thres = 0.0001;

err = err_Thres * ones(S_num,N);

max_err = err_Thres;

x_controlled_N = zeros(S_num, L,N);
for n = 1:N
    x_controlled_N(:,1,n)= P_default_N(S_num, :,N); % initial state
end

max_err_N = err_Thres * ones(N,1);



%% One-Stage Consensus
l = 1
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i, j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 2
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i, j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 3
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 4
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i,j, n) * x_controlled_N(j, l ,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 5
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 6
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 7
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 8
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 9
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 10
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State( i, j,n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

%% One-Stage Consensus
l = 11
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n) * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 12
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 13
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 14
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 15
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 16
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 17
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 18
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end
%% One-Stage Consensus
l = 19
Vec_State_temp = Vec_State(:,l,n);
for n = 1:N
    Matrix_State(:,:,n) = reshape(Vec_State_temp,[S_num,S_num]);
end
max_err = err_Thres;
while max_err >= err_Thres
    for n = 1:N
        
        Matrix_State_temp = zeros(S_num,S_num);
        
        for m = 1:N
            Matrix_State_temp(:,:) = G_mat(m,n) .* Matrix_State(:,:,m) + Matrix_State_temp(:,:) ;
        end
        epsilon = 0.000000000001;
        for i = 1:S_num
            for j = 1:S_num
                if Matrix_State_temp(i, j) == 0
                    Matrix_State_temp(i, j) = epsilon; % a small probability
                end
            end
        end
        for i = 1:S_num
            err(i,n) = sum(abs(Matrix_State(i,:,n)- Matrix_State_temp(i,:)));
            if err(i,n) <= err_Thres
            else
                for j = 1:S_num
                    D(i,j,n) = gamma(n)*(log(Matrix_State_temp(i,j)) - log(P_default_N(i,j,n))+1)+V(j,l+1,n);
                end
                Matrix_State_new(i,:,n) = Simplex_proj_vec(Matrix_State_temp(i,:) - 1/t*D(i,:,n),S_num);
                Matrix_State(i,:,n) = Matrix_State_new(i,:,n);
            end
        end
        max_err_N(n) = max(err(:,n));
        max_err = max(max_err_N);
        for i = 1:S_num
            temp3 = 0;
            for j = 1:S_num
                temp3 = temp3 + Matrix_State(i,j, n)  * x_controlled_N(j, l,n);
            end
            x_controlled_N(i, l+1, n) = temp3;
        end
    end
    
    t = t + 1; % next communication
end

for l = 1:L
    for n = 1:N
        elec_use_controlled_N(l,n) = elec_use * x_controlled_N(:,l,n);
    end
end

save('elec_use_controlled_N')
