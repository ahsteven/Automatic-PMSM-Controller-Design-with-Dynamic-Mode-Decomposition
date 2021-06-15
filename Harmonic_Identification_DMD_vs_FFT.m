clear all;

% load uncompensated data from Surface Mount PMSM
SPM_no_comp = load('./data/SPM_1000Rad_uncomp.mat');
SPM_no_comp = SPM_no_comp.data;

% load DMD compensated data from Surface Mount PMSM
SPM_DMD_comp = load('./data/SPM_1000Rad_DMD_comp.mat');
SPM_DMD_comp = SPM_DMD_comp.data;

% load ASHE compensation data from Surface Mount PMSM
% The ASHE algorithm is used for benchmark comparison
% Blasko, Vladimir. "A novel method for selective harmonic elimination 
% in power electronic equipment." IEEE transactions on Power 
% Electronics 22, no. 1 (2007): 223-228.
SPM_ASHE_comp = load('./data/SPM_1000Rad_ASHE_comp.mat');
SPM_ASHE_comp = SPM_ASHE_comp.data

% load uncompensated data from Interior Mount PMSM
IPM_no_comp = load('./data/IPM_200Rad_299delays_uncomp.mat');
IPM_no_comp = IPM_no_comp.data;

% load DMD compensated data from Interior Mount PMSM
IPM_DMD_comp = load('./data/IPM_200Rad_299delays_comp.mat');
IPM_DMD_comp = IPM_DMD_comp.data;

%% Calculate SPM Phase a Total Harmonic Distortion
% SPM 
SPM_no_comp_THD = thd(SPM_no_comp.ia1);% no compensation
SPM_DMD_comp_THD = thd(SPM_DMD_comp.ia1);% compensation with DMD 
SPM_ASHE_comp_THD = thd(SPM_ASHE_comp.ia1);% compensation with ASHE

% IPM
IPM_no_comp_THD = thd(IPM_no_comp.ia1);% no compensation
IPM_DMD_comp_THD = thd(IPM_DMD_comp.ia1);% no compensation

%% Calculate most dominant mode frequencies using DMD
%     DMD code from Scott Dawson see:
%     "Characterizing and correcting for the effect of sensor noise in the
%     dynamic mode decomposition"

% Setup DMD parameters
method = 'fb'; % forward backward DMD
ndelays = 299;
r = 8;

% for iq axis subtract of reference to view only the harmonic modes
% for id axis the reference is zero
% SPM no compensation
[ iq_eigs_SPM_no_comp, iq_modes_SPM_no_comp,~,iq_S_SPM_no_comp, iq_x0_SPM_no_comp] = dmd(SPM_no_comp.iq1-SPM_no_comp.iq1_ref,method,ndelays,r);
[ id_eigs_SPM_no_comp, id_modes_SPM_no_comp,~,id_S_SPM_no_comp, id_x0_SPM_no_comp] = dmd(SPM_no_comp.id1,method,ndelays,r);

% SPM compensation with DMD
[ iq_eigs_SPM_DMD_comp, iq_modes_SPM_DMD_comp,~,iq_S_SPM_DMD_comp, iq_x0_SPM_DMD_comp] = dmd(SPM_DMD_comp.iq1-SPM_DMD_comp.iq1_ref,method,ndelays,r);
[ id_eigs_SPM_DMD_comp, id_modes_SPM_DMD_comp,~,id_S_SPM_DMD_comp, id_x0_SPM_DMD_comp] = dmd(SPM_DMD_comp.id1,method,ndelays,r);

% SPM compensation with ASHE
[ iq_eigs_SPM_ASHE_comp, iq_modes_ASHE_comp,~,iq_S_SPM_ASHE_comp, iq_x0_ASHE_comp] = dmd(SPM_ASHE_comp.iq1-SPM_ASHE_comp.iq1_ref,method,ndelays,r);
[ id_eigs_SPM_ASHE_comp, id_modes_ASHE_comp,~,id_S_SPM_ASHE_comp, id_x0_ASHE_comp] = dmd(SPM_ASHE_comp.id1,method,ndelays,r);

% IPM no compensation
[ iq_eigs_IPM_no_comp, iq_modes_IPM_no_comp,~,iq_S_IPM_no_comp, iq_x0_IPM_no_comp] = dmd(IPM_no_comp.iq1-IPM_no_comp.iq1_ref,method,ndelays,r);
[ id_eigs_IPM_no_comp, id_modes_IPM_no_comp,~,id_S_IPM_no_comp, id_x0_IPM_no_comp] = dmd(IPM_no_comp.id1,method,ndelays,r);

% IPM DMD compensation
[ iq_eigs_IPM_DMD_comp, iq_modes_IPM_DMD_comp,~,iq_S_IPM_DMD_comp, iq_x0_IPM_DMD_comp] = dmd(IPM_DMD_comp.iq1-IPM_DMD_comp.iq1_ref,method,ndelays,r);
[ id_eigs_IPM_DMD_comp, id_modes_IPM_DMD_comp,~,id_S_IPM_DMD_comp, id_x0_IPM_DMD_comp] = dmd(IPM_DMD_comp.id1,method,ndelays,r);

%% View singular values to see dominant modes
% SPM DMD iq axis
figure()
stem(diag(iq_S_SPM_no_comp),'LineWidth',2); hold on
stem(diag(iq_S_SPM_DMD_comp),'LineWidth',2);
legend('Uncompensated','Compensated')
xlabel('SVD modes','Fontsize', 18)
ylabel('singular values','Fontsize', 18)
title('SPM Iq Singular Values Uncompensated Vs. DMD Compensated')  

% SPM DMD id axis
figure()
stem(diag(id_S_SPM_no_comp),'LineWidth',2); hold on
stem(diag(id_S_SPM_DMD_comp),'LineWidth',2);
legend('Uncompensated','Compensated')
xlabel('SVD modes','Fontsize', 18)
ylabel('singular values','Fontsize', 18)
title('SPM Id Singular Values Uncompensated Vs. DMD Compensated') 

% SPM ASHE iq axis
figure()
stem(diag(iq_S_SPM_no_comp),'LineWidth',2); hold on
stem(diag(iq_S_SPM_ASHE_comp),'LineWidth',2);
legend('Uncompensated','Compensated')
xlabel('SVD modes','Fontsize', 18)
ylabel('singular values','Fontsize', 18)
title('SPM Iq Singular Values Uncompensated Vs. ASHE Compensated') 

% SPM ASHE id axis
figure()
stem(diag(id_S_SPM_no_comp),'LineWidth',2); hold on
stem(diag(id_S_SPM_ASHE_comp),'LineWidth',2);
legend('Uncompensated','Compensated')
xlabel('SVD modes','Fontsize', 18)
ylabel('singular values','Fontsize', 18)
title('SPM Id Singular Values Uncompensated Vs. ASHE Compensated')

% View singular values for IPM
% % iq axis was not compensated for IPM
% % IPM DMD iq axis
% figure()
% stem(diag(iq_S_IPM_no_comp),'LineWidth',2); hold on
% stem(diag(iq_S_IPM_DMD_comp),'LineWidth',2);
% legend('Uncompensated','Compensated')
% xlabel('SVD modes','Fontsize', 18)
% ylabel('singular values','Fontsize', 18)
% title('IPM Iq Singular Values Uncompensated Vs. DMD Compensated') 

% IPM DMD id axis
figure()
stem(diag(id_S_IPM_no_comp),'LineWidth',2); hold on
stem(diag(id_S_IPM_DMD_comp),'LineWidth',2);
legend('Uncompensated','Compensated')
xlabel('SVD modes','Fontsize', 18)
ylabel('singular values','Fontsize', 18)
title('IPM Id Singular Values Uncompensated Vs. DMD Compensated') 

%% Stem plots of FFT and DMD frequencies
deltaT = SPM_no_comp.time(end)-SPM_no_comp.time(end-1); % switching frequecy 20kHz
% SPM uncompensated iq
[freqs_DMD_iq_SPM_no_comp, mags_DMD_iq_SPM_no_comp] = calc_DMD_freqs_mags(iq_eigs_SPM_no_comp, iq_modes_SPM_no_comp, iq_x0_SPM_no_comp, deltaT);
[freqs_fft_iq_SPM_no_comp, mags_fft_iq_SPM_no_comp] = calc_fft_freqs_mags(SPM_no_comp.iq1-SPM_no_comp.iq1_ref, deltaT);
% SPM uncompensated id
[freqs_DMD_id_SPM_no_comp, mags_DMD_id_SPM_no_comp] = calc_DMD_freqs_mags(id_eigs_SPM_no_comp, id_modes_SPM_no_comp, id_x0_SPM_no_comp, deltaT);
[freqs_fft_id_SPM_no_comp, mags_fft_id_SPM_no_comp] = calc_fft_freqs_mags(SPM_no_comp.id1, deltaT);
% SPM compensated with DMD iq
[freqs_DMD_iq_SPM_DMD_comp, mags_DMD_iq_SPM_DMD_comp] = calc_DMD_freqs_mags(iq_eigs_SPM_DMD_comp, iq_modes_SPM_DMD_comp, iq_x0_SPM_DMD_comp, deltaT);
[freqs_fft_iq_SPM_DMD_comp, mags_fft_iq_SPM_DMD_comp] = calc_fft_freqs_mags(SPM_DMD_comp.iq1-SPM_DMD_comp.iq1_ref, deltaT);
% SPM compensated with DMD id
[freqs_DMD_id_SPM_DMD_comp, mags_DMD_id_SPM_DMD_comp] = calc_DMD_freqs_mags(id_eigs_SPM_DMD_comp, id_modes_SPM_DMD_comp, id_x0_SPM_DMD_comp, deltaT);
[freqs_fft_id_SPM_DMD_comp, mags_fft_id_SPM_DMD_comp] = calc_fft_freqs_mags(SPM_DMD_comp.id1, deltaT);
% SPM compensated with ASHE iq
[freqs_DMD_iq_SPM_ASHE_comp, mags_DMD_iq_SPM_ASHE_comp] = calc_DMD_freqs_mags(iq_eigs_SPM_ASHE_comp, iq_modes_SPM_ASHE_comp, iq_x0_SPM_ASHE_comp, deltaT);
[freqs_fft_iq_SPM_ASHE_comp, mags_fft_iq_SPM_ASHE_comp] = calc_fft_freqs_mags(SPM_ASHE_comp.iq1-SPM_ASHE_comp.iq1_ref, deltaT);
% SPM compensated with ASHE id
[freqs_DMD_id_SPM_ASHE_comp, mags_DMD_id_SPM_ASHE_comp] = calc_DMD_freqs_mags(id_eigs_SPM_ASHE_comp, id_modes_SPM_ASHE_comp, id_x0_SPM_ASHE_comp, deltaT);
[freqs_fft_id_SPM_ASHE_comp, mags_fft_id_SPM_ASHE_comp] = calc_fft_freqs_mags(SPM_ASHE_comp.id1, deltaT);

figure()
stem(freqs_fft_iq_SPM_no_comp, mags_fft_iq_SPM_no_comp); hold on 
stem(abs(freqs_DMD_iq_SPM_no_comp), mags_DMD_iq_SPM_no_comp, 'linewidth', 2);
stem(freqs_fft_iq_SPM_DMD_comp, mags_fft_iq_SPM_DMD_comp); hold on 
% uses DMD for both the compensation and the post processing analysis
stem(abs(freqs_DMD_iq_SPM_DMD_comp), mags_DMD_iq_SPM_DMD_comp, 'linewidth', 2);
legend('FFT uncompensated','DMD uncompensated', 'FFT DMD compensated', 'DMD DMD compensated','fontsize', 18)
tit = ['Iq Current Uncompensated vs. DMD Compensated'];
title(tit)
ylabel('harmonic amplitude (A)','fontsize',18)
xlabel('frequency (rad/s)','fontsize',18) 

figure()
stem(freqs_fft_id_SPM_no_comp, mags_fft_id_SPM_no_comp); hold on 
stem(abs(freqs_DMD_id_SPM_no_comp), mags_DMD_id_SPM_no_comp, 'linewidth', 2);
stem(freqs_fft_id_SPM_DMD_comp, mags_fft_id_SPM_DMD_comp); hold on 
% uses DMD for both the compensation and the post processing analysis
stem(abs(freqs_DMD_id_SPM_DMD_comp), mags_DMD_id_SPM_DMD_comp, 'linewidth', 2);
legend('FFT uncompensated','DMD uncompensated', 'FFT DMD compensated', 'DMD DMD compensated','fontsize', 18)
tit = ['Id Current Uncompensated vs. DMD Compensated'];
title(tit)
ylabel('harmonic amplitude (A)','fontsize',18)
xlabel('frequency (rad/s)','fontsize',18) 

figure()
stem(freqs_fft_iq_SPM_no_comp, mags_fft_iq_SPM_no_comp); hold on 
stem(abs(freqs_DMD_iq_SPM_no_comp), mags_DMD_iq_SPM_no_comp, 'linewidth', 2);
stem(freqs_fft_iq_SPM_ASHE_comp, mags_fft_iq_SPM_ASHE_comp); hold on 
% uses ASHE for the compensation but DMD for post processing analysis
stem(abs(freqs_DMD_iq_SPM_ASHE_comp), mags_DMD_iq_SPM_ASHE_comp, 'linewidth', 2);
legend('FFT uncompensated','DMD uncompensated', 'FFT ASHE compensated', 'DMD ASHE compensated','fontsize', 18)
tit = ['Iq Current Uncompensated vs. ASHE Compensated'];
title(tit)
ylabel('harmonic amplitude (A)','fontsize',18)
xlabel('frequency (rad/s)','fontsize',18) 

figure()
stem(freqs_fft_id_SPM_no_comp, mags_fft_id_SPM_no_comp); hold on 
stem(abs(freqs_DMD_id_SPM_no_comp), mags_DMD_id_SPM_no_comp, 'linewidth', 2);
stem(freqs_fft_id_SPM_ASHE_comp, mags_fft_id_SPM_ASHE_comp); hold on 
% uses ASHE for the compensation but DMD for post processing analysis
stem(abs(freqs_DMD_id_SPM_ASHE_comp), mags_DMD_id_SPM_ASHE_comp, 'linewidth', 2);
legend('FFT uncompensated','DMD uncompensated', 'FFT ASHE compensated', 'DMD ASHE compensated','fontsize', 18)
tit = ['Id Current Uncompensated vs. ASHE Compensated'];
title(tit)
ylabel('harmonic amplitude (A)','fontsize',18)
xlabel('frequency (rad/s)','fontsize',18) 

function[freqs, mags] = calc_DMD_freqs_mags(DMD_eigs, DMD_vecs, x0, deltaT)
    cont_evals = log(DMD_eigs)/deltaT;
    % sort modes by magnitude
    b = pinv(DMD_vecs)*x0;
    [~, idx]=sort(abs(DMD_vecs(1,:)'.*b),'descend');

    sorted_evals = cont_evals(idx);
    sorted_vecs = DMD_vecs(:,idx);

    dmd_mags = abs(sorted_vecs(end,:)'.*pinv(sorted_vecs)*x0);
    dmd_freqs = abs(round(imag(sorted_evals)));

    % remove duplicate modes -- complex conjugage 
    freqs = dmd_freqs(1:2:end);
    mags = 2*dmd_mags(1:2:end);

end

function [freqs, mags] = calc_fft_freqs_mags(i1, deltaT)
    N = length(i1);
    freqs = 2*pi*(1/deltaT).*(0:(N/2))/N;
    siq_fft = fft(i1);
    mags = abs(siq_fft(1:length(freqs))/N);
    mags(2:end-1) = mags(2:end-1)*2;
end
