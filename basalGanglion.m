%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Logs       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Log 14th Jan 2020
%Define populations numbers

%%Log 16th Jan 2020
%Define connections weights

%%Log 6th Feb 2020
%Define mask arrays and separate model parameter arrays

%%Log 13th Feb 2020
%Define global parameter arrays referring LGN Matlab code and experiment.

%%Log 18th Feb 2020
%Define final probabilistic connections and color coding in plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Runtime
% TotalDuration = 500;
% TimeInt = 0.1;
TotalDataPoints = 60000;

%% Number of neurons
N_SNr = 27; 
N_GPe = 46;
N_STN = 14;
N_Str_FSI = 86;
N_Str_MSN_D1 = 50;
N_Str_MSN_D2 = 50;
N = N_Str_FSI + N_Str_MSN_D1 + N_Str_MSN_D2 + N_STN + N_GPe + N_SNr;

%% Model Parameters for STR-D1
strd1_a=0.02 * ones(N_Str_MSN_D1,1);
strd1_b=0.2 * ones(N_Str_MSN_D1,1);
strd1_c=-65.0 * ones(N_Str_MSN_D1,1);
strd1_d=8.0 * ones(N_Str_MSN_D1,1);
strd1_v_init = -80.0 * ones(N_Str_MSN_D1,1);
strd1_u_init = strd1_b .* strd1_v_init .* ones(N_Str_MSN_D1,1);

%% Model Parameters for STR-D2
strd2_a=0.02 * ones(N_Str_MSN_D2,1);
strd2_b=0.2 * ones(N_Str_MSN_D2,1);
strd2_c=-65.0 * ones(N_Str_MSN_D2,1);
strd2_d=8.0 * ones(N_Str_MSN_D2,1);
strd2_v_init = -80.0 * ones(N_Str_MSN_D2,1);
strd2_u_init = strd2_b .* strd2_v_init .* ones(N_Str_MSN_D2,1);

current_bias_str = -30.0;

%% Model Parameters for STR-FSI
fsi_a=0.1 * ones(N_Str_FSI,1);
fsi_b=0.2 * ones(N_Str_FSI,1);
fsi_c=-65.0 * ones(N_Str_FSI,1);
fsi_d=8.0 * ones(N_Str_FSI,1);
fsi_v_init = -70.0 * ones(N_Str_FSI,1);
fsi_u_init = fsi_b .* fsi_v_init .* ones(N_Str_FSI,1);

current_bias_fsi = -10.0;

%% Model Parameters for GPe
gpe_a=0.005 * ones(N_GPe,1);
gpe_b=0.585 * ones(N_GPe,1);
gpe_c=-65.0 * ones(N_GPe,1);
gpe_d=4.0 * ones(N_GPe,1);
gpe_v_init = -70.0 * ones(N_GPe,1);
gpe_u_init = gpe_b .* gpe_v_init .* ones(N_GPe,1);

current_bias_gpe = 2.0;

%% Model Parameters for SNR
snr_a = 0.005 * ones(N_SNr,1);
snr_b = 0.32 * ones(N_SNr,1);
snr_c = -65.0 * ones(N_SNr,1);
snr_d = 2.0 * ones(N_SNr,1);
snr_v_init = -70.0 * ones(N_SNr,1);
snr_u_init = snr_b .* snr_v_init .* ones(N_SNr,1);

current_bias_snr = 5.0;

%% Model Parameters for STN
stn_a=0.005 * ones(N_STN,1);
stn_b=0.265 * ones(N_STN,1);
stn_c=-65.0 * ones(N_STN,1);
stn_d=2.0 * ones(N_STN,1);
stn_v_init = -60.0 * ones(N_STN,1);
stn_u_init = stn_b .* stn_v_init .* ones(N_STN,1);

current_bias_stn = 5.0;

%% Conductance and Delay Params
tau_ampa = 6;
E_ampa = 0.0;
g_ampa = 0.5;

tau_gabba = 6.0;
E_gabba = -80.0;
g_gaba = 0.25;

mod_ampa_d2 = 0.2;
phi_msn_D = 2.75;
phi_fsi_D = 3.75;
phi_stn_D = -2;

%% Interconnections
mod_gaba = 0.073;
g_strd12snr = g_gaba * (1 + mod_gaba * phi_msn_D);
g_strd22snr = g_gaba * (1 - mod_gaba * phi_msn_D);
g_str2str = (1.0/2.55) * g_gaba;
g_gaba_gpe = (1.0/1.75) * g_gaba;
g_gpe2stn = g_gaba_gpe;
g_gpe2gpe = g_gaba_gpe;
g_gpe2snr = g_gaba_gpe;
g_gpe2fsi = g_gaba_gpe;
g_gaba_snr = (1/1.75) * g_gaba;
g_snr2snr = g_gaba_snr;
g_stn2str = (g_ampa * (1 - (mod_ampa_d2 * phi_stn_D))) / 6.0;
g_stn2gpe = (g_ampa * (1 - (mod_ampa_d2 * phi_stn_D))) / 6.0;

% a = [0.02*ones(N_Str_MSN_D1,1); 0.02*ones(N_Str_MSN_D2,1); 0.02*ones(N_Str_FSI, 1); 0.02*oneS(N_STN,1); 0.02*ones(N_SNr,1); 0.02*ones(N_GPe,1), 0.02*ones(N_)]
%STRD1, STRD2, STR_FSI, STN, GPE, SNR
a = [strd1_a;strd2_a;fsi_a;stn_a;gpe_a;snr_a];
b = [strd1_b;strd2_b;fsi_b;stn_b;gpe_b;snr_b];
c = [strd1_c;strd2_c;fsi_c;stn_c;gpe_c;snr_c];
d = [strd1_d;strd2_d;fsi_d;stn_d;gpe_d;snr_d];
v = [strd1_v_init;strd2_v_init;fsi_v_init;stn_v_init;gpe_v_init;snr_v_init];
u = [strd1_u_init;strd2_u_init;fsi_u_init;stn_u_init;gpe_u_init;snr_u_init];

% zeros array for storing connections
S = [zeros(N, N_Str_MSN_D1), zeros(N, N_Str_MSN_D2), zeros(N, N_Str_FSI), zeros(N, N_STN), zeros(N, N_GPe), zeros(N, N_SNr)];

% connection weights
S(1:N_Str_MSN_D1, 1:N_Str_MSN_D1) = g_str2str;% D1 to D2
S(1:N_Str_MSN_D1, 1+N_Str_MSN_D1:N_Str_MSN_D2) = g_str2str;% D1 to D2
S(1+N_Str_MSN_D1:N_Str_MSN_D2, 1:N_Str_MSN_D1) = g_str2str;% D2 to D1
S(1+N_Str_MSN_D1:N_Str_MSN_D2, 1+N_Str_MSN_D1:N_Str_MSN_D2) = g_str2str;% D2 to D2
S(1:N_Str_MSN_D1, 1+N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI+N_STN+N_GPe:N) = 0.3; %STR_D1 to SNr
S(1:N_Str_MSN_D1, 1+N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI+N_STN:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI+N_STN+N_GPe) = 0.2; %STR_D1 to GPe
S(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1:N_Str_MSN_D1) = 0.0982; %FSI to D1
S(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2) = 0.0982; % FSI to D2
S(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI) = 0.0982; % FSI to FSI
S(N-N_SNr+1:N, N-N_SNr+1:N) = 0.1429; %SNr to SNr
S(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe) = 0.1429; %GPe to STN
S(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr+1:N) = 0.1429; %GPe to SNr
S(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe+1:N-N_SNr) = 0.1429; %GPe to GPe
S(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe-N_STN-N_Str_FSI+1:N-N_SNr-N_STN-N_GPe) = 0.1429; %GPe to FSI
S(N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe, N-N_SNr-N_GPe+1:N-N_SNr) = 0.05; %STN to GPe
S(N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe, N-N_SNr+1:N) = 0.05; %STN to SNr

firings=[];             % spike timings

%% Defining masks
masks = [ones(N, N_Str_MSN_D1), ones(N, N_Str_MSN_D2), ones(N, N_Str_FSI), ones(N, N_STN), ones(N, N_GPe), ones(N, N_SNr)];

%Probabilistic connections defination. These masks are multiplied with the synaptic weights before updating neuron parameters
masks = update_masks(masks, N, N_Str_MSN_D1, N_Str_MSN_D2, N_Str_FSI, N_STN, N_GPe, N_SNr);

for t=1:TotalDataPoints            % simulation of 1000 ms
  I=[-30*randn(N_Str_MSN_D1,1); -30*randn(N_Str_MSN_D2,1); -10*randn(N_Str_FSI,1); 5*randn(N_STN,1); 2*randn(N_GPe,1); 5*randn(N_SNr, 1)]; % thalamic input
  fired=find(v(:,t)>=30);    % indices of spikes
  firings=[firings; t+0*fired,fired];
  v(fired,t)=c(fired);
  u(fired,t)=u(fired,t)+d(fired);

  temp_S = S .* masks; %Create temporary synaptic weight array after considering probabilistic connections.
  masks = update_masks(masks, N, N_Str_MSN_D1, N_Str_MSN_D2, N_Str_FSI, N_STN, N_GPe, N_SNr); %Create new masks for next loop
  I=I+sum(temp_S(:,fired),2); %Update I

  v(:,t+1)=(v(:,t)+0.05*(0.04*v(:,t).^2+5*v(:,t)+140-u(:,t)+I)); % step 0.5 ms
  v(:,t+1)=(v(:,t+1)+0.05*(0.04*v(:,t+1).^2+5*v(:,t+1)+140-u(:,t)+I)); % for numerical
  u(:,t+1)=u(:,t)+0.1*a.*(b.*v(:,t)-u(:,t));                 % stability
end;

%% Plotting;
y1 = firings(:,2);y2 = firings(:,2); y3 = firings(:,2); y4 = firings(:,2); y5 = firings(:,2); y6 = firings(:,2);
y1(firings(:,2)>51) = NaN;                                                                                      
y2(firings(:,2)>101 | firings(:,2)<50) = NaN;                                                                   
y3(firings(:,2)>187 | firings(:,2)<100) = NaN;                                                                  
y4(firings(:,2)>201 | firings(:,2)<186) = NaN;                                                                  
y5(firings(:,2)>247 | firings(:,2)<200) = NaN;                                                                  
y6(firings(:,2)>274 | firings(:,2)<246) = NaN;
yi = [y1, y2, y3, y4, y5, y6];linecolor = ['b' 'g' 'k' 'c' 'm' 'r'];
figure(1);
hold on;
for i=1:6
  % figure(1);
  plot(firings(:,1), yi(:,i), '.', 'Color', linecolor(i), 'MarkerSize',3);
end
legend({'STR_MSN_D1','STR_MSN_D2', 'STR_FSI', 'STN', 'GPe', 'SNr'}, 'Location', 'bestoutside');title("Rasterplot");
%plot(firings(:,1),firings(:,2),'.');title("Rasterplot");%annotation('textbox',[0 .5 .1 .2],'String',"1:" + N_Str_MSN_D1 + "STR-MSN D1\n" + N_Str_MSN_D1 + ":" + (N_Str_MSN_D1+N_Str_MSN_D2) + "STR-MSN D2\n",'EdgeColor','none');

figure;
subplot(6,1,1);plot(1:60001,mean(v(1:N_Str_MSN_D1,:)));title('Mean Neuron Potential of STR MSN D1');
subplot(6,1,2);plot(1:60001,mean(v(N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2,:)));title('Mean Neuron Potential of STR MSN D2');
subplot(6,1,3);plot(1:60001,mean(v(N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI,:)));title('Mean Neuron Potential of STR MSN FSI');
subplot(6,1,4);plot(1:60001,mean(v(N-N_SNr-N_GPe-N_STN:N-N_SNr-N_GPe,:)));title('Mean Neuron Potential of STN');
subplot(6,1,5);plot(1:60001,mean(v(N-N_SNr-N_GPe:N-N_SNr,:)));title('Mean Neuron Potential of GPe');
subplot(6,1,6);plot(1:60001,mean(v(N-N_SNr:N,:)));title('Mean Neuron Potential of SNr');

function masks = update_masks(masks, N, N_Str_MSN_D1, N_Str_MSN_D2, N_Str_FSI, N_STN, N_GPe, N_SNr)
  masks(1:N_Str_MSN_D1, 1:N_Str_MSN_D1) = rand(N_Str_MSN_D1, N_Str_MSN_D1) < 0.1; %Connection from D1 to D2 with probability 0.1
  masks(1:N_Str_MSN_D1, 1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2) = rand(N_Str_MSN_D1, N_Str_MSN_D2) < 0.1; %Connection from D1 to D2 with probability 0.1
  masks(1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2, 1:N_Str_MSN_D1) = rand(N_Str_MSN_D1, N_Str_MSN_D2) < 0.1; %Connection from D2 to D1 with probability 0.1
  masks(1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2, 1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2) = rand(N_Str_MSN_D1, N_Str_MSN_D2) < 0.1; %Connection from D2 to D2 with probability 0.1
  masks(1:N_Str_MSN_D1, N-N_SNr+1:N) = rand(N_Str_MSN_D1, N_SNr) < 0.15; %D1 to SNr
  masks(1:N_Str_MSN_D1, N-N_SNr-N_GPe+1:N-N_SNr) = rand(N_Str_MSN_D1, N_GPe) < 0.15; %D1 to GPe
  masks(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1:N_Str_MSN_D1) = rand(N_Str_FSI, N_Str_MSN_D1) < 0.1; %FSI to D1
  masks(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1+N_Str_MSN_D1:N_Str_MSN_D1+N_Str_MSN_D2) = rand(N_Str_FSI, N_Str_MSN_D2) < 0.1; % FSI to D2
  masks(1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI, 1+N_Str_MSN_D1+N_Str_MSN_D2:N_Str_MSN_D1+N_Str_MSN_D2+N_Str_FSI) = rand(N_Str_FSI, N_Str_FSI) < 0.1;
  masks(N-N_SNr+1:N, N-N_SNr+1:N) = rand(N_SNr, N_SNr) < 0.25; %SNr to SNr
  masks(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe) = rand(N_GPe, N_STN) < 0.25; %GPe to STN
  masks(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr+1:N) = rand(N_GPe, N_SNr) < 0.25; %GPe to SNr
  masks(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe+1:N-N_SNr) = rand(N_GPe, N_GPe) < 0.25; %GPe to GPe
  masks(N-N_SNr-N_GPe+1:N-N_SNr, N-N_SNr-N_GPe-N_STN-N_Str_FSI+1:N-N_SNr-N_STN-N_GPe) = rand(N_GPe, N_Str_FSI) < 0.05; %GPe to FSI
  masks(N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe, N-N_SNr-N_GPe+1:N-N_SNr) = rand(N_STN, N_GPe) < 0.5;%STN to GPe
  masks(N-N_SNr-N_GPe-N_STN+1:N-N_SNr-N_GPe, N-N_SNr+1:N) = rand(N_STN, N_SNr) < 0.5; %STN to SNr
end