clear all; close all

%% This code demmonstrates the ASHE algorithm from the paper:
% Blasko, Vladimir. "A novel method for selective harmonic elimination 
% in power electronic equipment." IEEE transactions on Power 
% Electronics 22, no. 1 (2007): 223-228.

t = [0:.001:.8];
ws = 1;
wc = 1;
w0 = 50;
phase0 = .4;
noise = 0*rand(size(t));
signal = 1.3*cos(w0.*t+phase0)+noise;

x_cos = cos(w0.*t); % generate sine and cosine waves at known harmonic frequency
x_sin = sin(w0.*t);
error = 0;
mu = .01; % set learning rate for gradient descent
wk_cos = 0;
wk_sin = 0;
cos_weights = [];
sin_weights = [];
cos_adj = [];
sin_adj = [];
pred_signal = [];
%%
for ii = 1:length(t);
    grad_e = error*2*mu; % calculate gradient of the error
    wk_cos = grad_e*x_cos(ii) + wk_cos;% update weights via gradient descent
    wk_sin = grad_e*x_sin(ii) + wk_sin;
    cos_weights = [cos_weights wk_cos];
    sin_weights = [sin_weights wk_sin];
    cos_adj =[cos_adj x_cos(ii)*wk_cos];
    sin_adj =[sin_adj x_sin(ii)*wk_sin];
    comb =cos_adj(end) + sin_adj(end);% predicted harmonic component
    pred_signal = [pred_signal comb];
    error = signal(ii) - comb;
    
end


figure()
plot(t, signal,'r','Linewidth',2);hold on
plot(t, pred_signal,'c','Linewidth',2); hold on
plot(t,cos_adj,'b','Linewidth',2)
plot(t,cos_weights,'Color','b','LineStyle','--','Linewidth',2)
plot(t,sin_adj,'g','Linewidth',2)
plot(t,sin_weights,'Color','g','LineStyle','--','Linewidth',2)
legend('signal','predicted','adj cos','weight cos','adj sin','weight sin','FontSize',14,'NumColumns',3,'Location','Southwest')
title('ASHE')
axis([0,t(end),-2,1.5])

