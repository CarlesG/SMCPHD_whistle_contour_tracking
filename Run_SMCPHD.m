%% ~~~~~~~~~~~~ SMC-PHD for dolphin whistle tracking ~~~~~~~~~~~~~~~~~~~~~
clear

%This example runs the version of the SMC-PHD filter that uses the Radial
%Basis Function (RBF) network to draw persistent particles in the prediction step

%For details see Gruden, P. and White, P. (2020). Automated extraction of dolphin whistles - 
%a Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)
%approach; Journal of the Acoustical Society of America.

%% ////////////  LOAD AUDIO DATA                           //////////////
clear all, close all
[file, path] = uigetfile('*.wav', 'Selected the file to analyse...');
[x,fs] = audioread([path file]);
disp(['Working on the file: ' file])
x = x(:,1); %select first channel

% load('Measurement_and_GT_demo.mat');
slide_incr = 1024;
dt = 0.0053; %time increment in s (default)
win_width = 2048;
win_width_s = 2048/fs;
peak_thr = 6;
Nvalid = 292;
freqrange = [5000, 35000];

% GT = ground truth data
% Zset = measurements (spectral peaks)
% dt = time step between consecutive windows (s)

%% //////////// INITIALIZE SMC-PHD PARAMETERS and MODELS /////////////
% %--------------- Parameters for SMCPHD ------------------------------------
% Table 1 in Gruden & White (2020):
parameters.pdet = 0.99; %probability of detection 
parameters.psurv = 0.994;%probability of survival 
parameters.Mp=50; %number of particles per persistent target
parameters.Nb=50; %number of particles per newborn target
parameters.nClutter = 10; % number of clutter points per time step
parameters.wth=0.0005; % threshold for state estimation (\eta) 

% Load GMM distribution for drawing the chirp - used in target birth step
% Eq. (9) in Gruden & White (2020):
load('gmm_chirp'); 
parameters.gmm_all=gmm_all;

% %-------------- Models for SMCPHD---------------:
models.H = [1,0]; % measurement matrix
models.R = round((fs/win_width)^2/12); % measurement noise variance for Eq. (5) in Gruden & White (2020)
load('birthpdf.mat'); %start frequency distribution for drawing adaptive 
% weight in Eq. (10) in Gruden & White (2020)
models.birthpdf=birthpdf;
models.dt=dt;
models.F = [1, dt; 0, 1]; % state transition matrix (used to resolve state
%estimate label conflicts in  track_labels.m)

%--------------- Parameters for RBF network--------------------------------
% Learned RBF network- centres, variances and weights (from training data)
% Eq. (7):
load('Net_60RBF.mat')
RBFnet.C=C; % (c_j)
RBFnet.w=w; % (w_j)
RBFnet.vari=vari; % (Q_j)
% Learned process noise (n_k) - Eq. (7):
models.Q = [39.2,0;0,7326]; 

% Settings for split the audio file in sections of sec_analysis duration
sec_analysis = 1;
N_analysis = round(sec_analysis * fs);
N_max = length(x);
idx = 1:N_analysis:N_max;
E_sp = [];

for i = 1:length(idx)
%     clf(1)
    disp(['Iteration: ' num2str(i) ' of ' num2str(length(idx))])
    if i == length(idx)
        y_part = x(idx(i):end);
    else
        y_part = x(idx(i) : idx(i) + N_analysis - 1);
    end
    slide_incr = round(dt*fs);
    numstps = ceil((length(y_part) - win_width) / slide_incr); %Number of overlapping windows
    %% //////////// MAKE MEASUREMENTS (Zsets) ///////////////
    disp('Making measurements')
    [Zset] = preprocess_getZset(win_width_s,dt,fs,y_part,freqrange,peak_thr);

    %% /////////// RUN SMC-PHD filter /////////////////
    disp('Running SMC-PHD filter')
    %Get state estimates and their labels:
    [Xk,XkTag] = SMCPHD_RBF_adaptivebirth(Zset,parameters,models,RBFnet);
    
    % Collect estimated states into tracks (whistle contours) based on their labels:
    disp('Tracking based on their labels')
    Track  = track_labels(XkTag,Xk,models);
    
    % Impose track length criteria to get detections that match specified criterion:
    tl=10;%target length criteria
    c=1;ind=[];
    for l=1:size(Track,2)
        if numel(Track(l).time)>=tl
            ind(c)=l;
            c=c+1;
        else
            continue
        end
    end
    DT_sp = Track(ind); % detected whistles
    % Aqui faltaría comprobar si DT está vacío. Sino, no se crearía el
    % evento de silbido
        
    %% ////////// PLOT RESULTS ///////////////
    figure(1),clf
     %plot measurements
    t = (0:1:(size(Zset,2)-1)) * models.dt + sec_analysis * (i-1); %time vector for plotting
    for k=1:size(Zset,2)
        if ~isempty(Zset{k})
           plot(t(k),Zset{k},'k.'),hold on
        end
    end
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    xlim([min(t) max(t)])
    title('Measurements')
    
%      subplot(212)
%     % Plot GT:
%     col= [0,0,0]+0.8;
%     for k=1:size(GT,2)
%         if GT(k).valid==1
%         plot(GT(k).time,GT(k).freq,'-','Color', col,'LineWidth',4),hold on
%         else
%         plot(GT(k).time,GT(k).freq,':','Color', col,'LineWidth',4),hold on    
%         end
%     end
    
    % Plot Detections:
    figure(1)
    for k=1:size(DT_sp,2)
        DT_sp(k).time = DT_sp(k).time + idx(i)/fs;
        plot(DT_sp(k).time,DT_sp(k).freq,'-','LineWidth',2)
    end
    disp(['Number of whistles detected: ' num2str(length(DT_sp))])
    ylabel('Frequency (kHz)','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
%     xlim([0,T(end)])
    title('Candidates spectrogram and SMC-PHD detections','Interpreter','latex')
    yt = freqrange(1):2000:freqrange(2);
    set(gca,'YTick',yt,'YTicklabel',yt*1e-3);
    set(gcf,'color','w');box on; grid;
    % Spectrogram representation in dB/Hz
    figure(3),spectrogram(y_part,2048,[],freqrange(1):100:freqrange(2),fs,'yaxis'); %remember that the units are dB/Hz 
    E_sp = [E_sp DT_sp];
    DT_sp = [];
end
close all