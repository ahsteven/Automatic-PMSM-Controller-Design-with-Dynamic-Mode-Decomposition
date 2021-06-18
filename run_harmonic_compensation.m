% This file runs the harmonic compensation simulation
% It shows that the algorithm is able to compensate
% 6, 12 and 18 th order harmonics satisfactorily 
% even with changing speed

clear all; close all;

deltaT = .00005;
Nsim = 20000;

params.rs = .7;
params.Ld_dc = .1447;
params.Lq_dc = .1631;
params.lmda_f_dc = .4055/(2*pi);
params.Poles = 14;
params.TL = 0;

params.Bm = .00011;
params.J = 1e-2;
params.spd_ctl = true;

lam_pm = params.lmda_f_dc;
Lq = params.Lq_dc;
Ld = params.Ld_dc;
Rs = params.rs;
Bm = params.Bm;
J = params.J;
P = params.Poles;

w_tuned = 1;
Acont = [-params.rs/params.Lq_dc, -params.Ld_dc*w_tuned/params.Lq_dc;...
    w_tuned*params.Lq_dc/params.Ld_dc, -params.rs/params.Ld_dc];
Bcont= [1/params.Lq_dc, 0, -params.lmda_f_dc/params.Lq_dc; 0, 1/params.Ld_dc, 0];
Ad = expm(Acont*deltaT);
Bd = (expm(Acont*deltaT) - eye(2))*pinv(Acont)*Bcont;
Ad_int = [Ad, zeros(2,2);-eye(2)*deltaT,eye(2)];
Bd_int = [Bd(:,1:2); zeros(2,2)];% first two columns only
Kd0 = [1722.92514708863,-0.149487775533096,...
    -243.708496078677,0;...
    0.158174537169950,1528.50661836498,...
    0,-216.219080907873];
B13 = -params.lmda_f_dc;

% determine DMD compensation vectors
amps = [1,1,1];
orders = [6,12,18];
start_comp = 50;
speeds = start_comp:1:1000;
comp_table0 = [];
comp_table1 = [];
nD = 299;
for ii = 1:length(speeds)
    s = speeds(ii)
    h_freqs = orders*s;
%     [vals, vecs, inv_vecs, x0] = get_comp_vecs_multi(amps,h_freqs,1,T_sw, nD, r);
    [vals, vecs, inv_vecs] = get_comp_vecs_algo_multi(h_freqs, deltaT, nD);
    comp_table1 = [comp_table1; real(vecs(1,:)*diag(vals)*inv_vecs)];
end

% save('tables_50to1000.mat','speeds','comp_table1');
% load('tables_50to1000.mat','speeds','comp_table1');
%%
comp_table1_q = comp_table1;
comp_table1_d = comp_table1;

% feed forward gain
H = pinv(-Ad^-1*Bd);
ss_q = H(1,1);
ss_d = H(2,2);

max_predq = .15;


%% Run without compensatin

run_comp = 0;
SimOut = sim('./sim_harm_comp/harmonic_comensation_simulation.slx');
iqd_no=SimOut.get('iqd');
w_e=SimOut.get('w_e');

run_comp = 1;
SimOut = sim('./sim_harm_comp/harmonic_comensation_simulation.slx');
iqd=SimOut.get('iqd');
q_comp = SimOut.get('q_comp');
d_comp = SimOut.get('d_comp');

figure()
plot(iqd_no(:,1),iqd_no(:,2)); hold on
plot(iqd(:,1),iqd(:,2));
legend('iq no comp','iq comp');
title('Iq compensation vs non compensated')

figure()
plot(iqd_no(:,1),iqd_no(:,3)); hold on
plot(iqd(:,1),iqd(:,3));
legend('id no comp','id comp');
title('Iq compensation vs non compensated')
% figure()
% plot(q_comp(:,1),q_comp(:,2));hold on
% plot(d_comp(:,1),d_comp(:,2));hold on

figure()
plot(w_e(:,1), w_e(:,2));
title('Electrical speed rad/s')

