% Decentrailized Ensemble Control
% Guanze Peng
% Janurary, 2022

%% System Setup
clear all
close all

global N
N = 100; % Number of units
global L
L = 20; % Number of stages for the period 9 a.m. to 3 p.m. 15 min intervals
global S_num
S_num = 15; % Number of states

P_default = load('MP_temp_mine.csv'); % Default matrix
epsilon = 0.000000000001;
for i = 1:S_num
    for j = 1:S_num
        if P_default(i, j) == 0
            P_default(i, j) = epsilon; % a small probability
        end
    end
end


% Electricity usage in MW
elec_use = 11.6:(34.5-11.6)/(S_num-1):34.5;
elec_use = elec_use/100; % convert to kW

% Electricity price 1 usd/kW
price = [42.79*ones(4,1);  43.61* ones(4,1); 46.38 *ones(4,1); 46.91*ones(4,1); 47.00*ones(4,1)]; % day ahead price 06/02/2018, 11:00 - 15: 00


%% Generate each unit's default matrix and discomfort cost

% Generate Discomfort Price
global gamma
gamma = 20 * ones(N,1);
Range_gamma = 4;
for n = 1:N
    r_gamma = floor(Range_gamma*rand(1))-Range_gamma; % generate random row vector
    gamma(n) = gamma(n) + r_gamma;
end

% Generate random default matrix
beta = 0.95; % parameter 1-beta is the level of how much the units are subject to randomness
global P_default_N
P_default_N = zeros(S_num,S_num,N);
for n = 1:N
    for i = 1:S_num
        r(i,:) = rand(1,S_num); % generate random row vector
        r_temp(i) = sum(r(i,:));  % Normalize so the sum is 1.
    end

    for i = 1:S_num
        for j = 1:S_num
            P_default_N(i,j,n) = P_default(i,j) * beta + (1-beta) * r(i, j)/r_temp(i);
            if P_default_N(i,j,n) == 0
                P_default_N(i,j,n) = epsilon; % a small probability
            end
        end
    end
end
Initial_state = zeros(S_num,1);
for n = 1:N
    Initial_state = Initial_state + P_default_N(:,1,n);
end
Initial_state = Initial_state./N;

% Cost from electricity consumption
global ele_cost
ele_cost  = zeros(S_num,L);
for l = 1:L
    for i = 1:S_num
        ele_cost(i, l) =  elec_use(i) * price(l);
    end
end


%% Generate communication matrix
G = ones(N,N)*(1/N);

global G_mat
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

%% Generate communication matrix for global consensus
G = ones(N,N)*(1/N);

global G_mat
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

%% Execute

main_centrailized
main_decentrailized
main_global_consensus

load('elec_use_controlled')
load('elec_use_controlled_N')
load('elec_use_controlled_global')


figure(1)


% Invidual controlled
for l = 1:L
    elec_use_max(l) = max(elec_use_controlled(l,:));
    elec_use_min(l) = min(elec_use_controlled(l,:));
end

time_day = 0:24;
stage_day = 1:1:19*60/15;
stage_dr = stage_day(21:20+L);

h(1) = plot(stage_dr, elec_use_max, '-.b', 'LineWidth',1)
hold on
h(2) = plot(stage_dr, elec_use_min, '-.b', 'LineWidth',1)
time2 = [stage_dr, fliplr(stage_dr)];
inBetween = [elec_use_min, fliplr(elec_use_max)];
h(3) = fill(time2, inBetween, [0.7 0.8 1]);

% x_default_day_before = Initial_state;
% % Default ensemble vector
% for l = 1:20
%     for i = 1:S_num
%         temp3 = 0;
%         for j = 1:S_num
%             temp3 = temp3 + P_default(j,i) * x_default_day_before(j, l);
%         end
%         x_default_day_before(i, l+1) = temp3;
%     end
% end

% Electricity consumption (default and controlled)
for l = 1:20
    elec_use_default_before_temp(l) = 0.06 * exp(0.08*(l-1));
end
Initial_use = elec_use * Initial_state;
% elec_use_default_before = elec_use_default_before_temp(end:-1:1);
elec_use_default_before = elec_use_default_before_temp;
elec_use_default_all_temp = cat(2, elec_use_default_before(1:end), elec_use_default_ave);



x_default_day_after(:,1) = x_default_ave(:,end);
% Default ensemble vector
for l = 1:24*60/15-10.5*4-L+2
    for i = 1:S_num
        temp3 = 0;
        for j = 1:S_num
            temp3 = temp3 + P_default(j, i) * x_default_day_after(j, l);
        end
        x_default_day_after(i, l+1) = temp3;
    end
end

% Electricity consumption (default and controlled)
for l = 1:36
    elec_use_default_after(l) = elec_use * x_default_day_after(:,l+1);
end

elec_use_default_all = cat(2, elec_use_default_all_temp, elec_use_default_after);


h(4) = plot(stage_day, elec_use_default_all, '-b', 'LineWidth',2)
h(5) = plot(stage_dr, elec_use_controlled_bar, '-.r', 'LineWidth',2)
h(6) = plot(stage_dr, elec_use_controlled_N(:,1), '-.g', 'LineWidth',2)
h(7) = plot(stage_dr, elec_use_controlled_global(:,1), '-.c', 'LineWidth',2)

legend(h([4 5 6 7 3]),'Default Consumption','Centrailized Controller','Decentrailized Controller with Local Consensus','Decentrailized Controller with Global Consensus','Range of Optimal Consumption without Consensus')


set(gca,'FontSize',15)
grid on
xlabel('Time (min)')
ylabel('Power (MW)')

set(gca,'FontSize',15)
grid on
xlabel('Time (Hour)')
ylabel('Power (MW)')

% x_hour = linspace(0, 24, 24*60/15)
% set(gca,'XTick',x_hour)

% for n = 1 : N
%     plot(time, elec_use_default(:,n), '-k')
%     hold on
%     plot(time, elec_use_controlled(:,n), '-.r')
% end
