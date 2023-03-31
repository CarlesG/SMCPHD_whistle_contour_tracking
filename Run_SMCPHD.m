%% ~~~~~~~~~~~~ SMC-PHD for dolphin whistle tracking ~~~~~~~~~~~~~~~~~~~~~
clear

%This example runs the version of the SMC-PHD filter that uses the Radial
%Basis Function (RBF) network to draw persistent particles in the prediction step

%For details see Gruden, P. and White, P. (2020). Automated extraction of dolphin whistles - 
%a Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)
%approach; Journal of the Acoustical Society of America.

%% ////////////  LOAD AUDIO DATA                           //////////////
clear all, close all
addpath('pyknogram_functions'); % All pyknogram functions
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


% Settings for the definition of the filter Bank of piknogram
BW=1000; % BW in Hz
BWoverlap=50; % BW in %
flow=freqrange(1);fhigh=min([freqrange(2) fs/2-BW/2]); % MINIMUM flow -> flow=round(BW/2)
Tl = sec_analysis;

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
    %% //////////// PIKNOGRAM IMPLEMENTATION ////////////////

    [FW,BW_est, ndraw] = pyknogram_freqdomain(y_part, fs, flow, fhigh, BW, BWoverlap, 10e-3);
    X=repmat(ndraw,1,size(FW,2)); 
    %  Get points with a Bandwidth lower than a given value (BW_thre)
    BW_thre=BW/4;
    BW_est(1,:)=BW; % To get rid of some problems in the BW estimation of the first and last
    BW_est(end,:)=BW; % To get rid of some problems in the BW estimation of the first and last
    idx_BW=find(BW_est<BW_thre);
    
    % ----------
    %  Point density based filter
    [P,idx_new]  = kernel_density( FW, ndraw, BW, BWoverlap );
    idx_combined = intersect(idx_BW,idx_new');   
    
    % Create a struct and store possible whistles
    Pyk=struct('time',cell(1,length(idx_combined)),'freq',[],'ampl',[],'label',[],'done',[]);
    kk=num2cell(X(idx_combined));[Pyk.time]=kk{:};
    kk=num2cell(FW(idx_combined));[Pyk.freq]=kk{:};
    
    % Represent candidates
%     figure(2);clf;
%     %sh=scatter([Pyk.time],[Pyk.freq],[],[Pyk.ampl],'filled');
%     sh=scatter([Pyk.time],[Pyk.freq],[],'b','filled');
%     sh.SizeData=15; 
%     title(['PRWE method']);
%     ylabel('Frequency [kHz]');xlabel('Time [sec.]');
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([0 Tl flow fhigh]);
%     set(gcf,'color','w');box on;grid;

    % Ordering candidates;
    tiempo = [Pyk.time];
    frecuencia = [Pyk.freq];
    [tiempo_ordenado, ind_ordenado] = sort(tiempo);
    frecuencia_ordenado = frecuencia(ind_ordenado);
    [tiempo_unico,ia,ic] = unique(tiempo_ordenado);
    Zset = cell([1 length(tiempo_unico)]);
    ic_unico = unique(ic);

    for ii = 1:length(Zset)
        Zset{ii} = sort(frecuencia_ordenado(ic == ic_unico(ii)));
    end
    
    % Represent ordered candidates
%     figure(2),clf
%     for j = 1:length(tiempo_unico)
%         plot(tiempo_unico(j),Zset{j},'b.'),hold on
%     end
%     title('Candidates Picknogram')
%     ylabel('Frequency [kHz]');xlabel('Time [sec.]');
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([0 Tl flow fhigh]);
%     set(gcf,'color','w');box on;grid;

    % Conformation of all complet Zset candidates
    Zset_all = cell(1, numstps);
    t=(0:(size(Zset_all,2))).*dt;  
    
    for z = 1:length(Zset)
        [~,ix] = min(abs(t - tiempo_unico(z)));
        Zset_all{ix} = Zset{z};
    end

    Zset_pic = Zset_all;
%     %count the number of elements different of empty in the cell array
    %count = length(Zset_all) - nnz(cellfun(@isempty,Zset_all))
    %% //////////// MAKE MEASUREMENTS (Zsets) ///////////////
    disp('Making measurements')
    [Zset_sp] = preprocess_getZset(win_width_s,dt,fs,y_part,freqrange,peak_thr);

    %% /////////// RUN SMC-PHD filter /////////////////
    disp('Running SMC-PHD filter')
    %Get state estimates and their labels:
    [Xk_pic,XkTag_pic] = SMCPHD_RBF_adaptivebirth(Zset_pic,parameters,models,RBFnet);
    [Xk_sp,XkTag_sp] = SMCPHD_RBF_adaptivebirth(Zset_sp,parameters,models,RBFnet);
    
    % Collect estimated states into tracks (whistle contours) based on their labels:
    disp('Tracking based on their labels')
    Track_pic  = track_labels(XkTag_pic,Xk_pic,models);
    Track_sp  = track_labels(XkTag_sp,Xk_sp,models);
    
    % Impose track length criteria to get detections that match specified criterion:
    tl_sp=10;%target length criteria
    c_sp=1;ind_sp=[];
    for l=1:size(Track_sp,2)
        if numel(Track_sp(l).time)>=tl_sp
            ind_sp(c_sp)=l;
            c_sp=c_sp+1;
        else
            continue
        end
    end

    tl_pic=3;%target length criteria
    c_pic=1;ind_pic=[];
    for l=1:size(Track_pic,2)
        if numel(Track_pic(l).time)>=tl_pic
            ind_pic(c_pic)=l;
            c_pic=c_pic+1;
        else
            continue
        end
    end

    if ~isempty(ind_sp)
        DT_sp = Track_sp(ind_sp); % detected whistles
    else
        DT_sp = [];
    end

    if ~isempty(ind_pic)
        DT_pic = Track_pic(ind_pic); % detected whistles
    else
        DT_pic = [];
    end       
    %% ////////// PLOT RESULTS DETECTIONS SPECTROGRAM ///////////////
    figure(1),clf
     %plot measurements
    t = (0:1:(size(Zset_sp,2)-1)) * models.dt + sec_analysis * (i-1); %time vector for plotting
    for k=1:size(Zset_sp,2)
        if ~isempty(Zset_sp{k})
           plot(t(k),Zset_sp{k},'k.'),hold on
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
    for k=1:size(DT_sp,2)
        DT_sp(k).time = DT_sp(k).time + idx(i)/fs;
        plot(DT_sp(k).time,DT_sp(k).freq,'-','LineWidth',2)
    end
    disp(['Number of whistles detected: ' num2str(length(DT_sp))])
    ylabel('Frequency (kHz)','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
%     xlim([0,T(end)])
    title('Candidates vs detections SPECTROGRAM','Interpreter','latex')
    yt = freqrange(1):2000:freqrange(2);
    set(gca,'YTick',yt,'YTicklabel',yt*1e-3);
    set(gcf,'color','w');box on; grid;

    figure(2),clf
     %plot measurements
    t = (0:1:(size(Zset_pic,2)-1)) * models.dt + sec_analysis * (i-1); %time vector for plotting
    for k=1:size(Zset_pic,2)
        if ~isempty(Zset_pic{k})
           plot(t(k),Zset_pic{k},'k.'),hold on
        end
    end
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    xlim([min(t) max(t)])
%     title('Measurements')
    
    for m=1:size(DT_pic,2)
	    DT_pic(m).time = (DT_pic(m).time * fs + idx(i))./fs;
        plot(DT_pic(m).time, DT_pic(m).freq,'LineWidth',1.5),hold on
    end

    title('Candidates vs detections PIKNOGRAM','Interpreter','latex')
    xlabel('time (s)','Interpreter','latex')
    ylabel('frequency (kHz)','Interpreter','latex')
    yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    axis([min(t) max(t) flow fhigh]);
    set(gcf,'color','w');box on;grid;


    % Spectrogram representation in dB/Hz
    figure(3),spectrogram(y_part,2048,[],freqrange(1):100:freqrange(2),fs,'yaxis'); %remember that the units are dB/Hz 
    E_sp = [E_sp DT_sp];
    DT_sp = [];
end
close all