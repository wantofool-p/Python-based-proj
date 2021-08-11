clear;
PSD_dBm=-174;
bandwidth=20e6;
d1_km=0.1;
d2_km=0.5;
N=6; %transmitter
K=4; %receiver
Pmax_dBm=30;

PSD_linear=10^((PSD_dBm-30)/10);
noise_pow=PSD_linear*bandwidth;
[H,w,n_complex]=hk_producer(d1_km,d2_km,N,K);

Pmax_linear=10^((Pmax_dBm-30)/10);
% W=zeros(N,K);
% for user=1:K
%     w=sqrt(1/2)*(randn(N,1)+randn(N,1)*1i);
%     w=w/norm(w);
%     w=w*sqrt(Pmax_linear/K);
%     W(:,user)=w;
% end    
[W,D]=ZFBF(H,Pmax_dBm);
SINR=getSINR(H,W,noise_pow)
[W1,D1]=RZFBF(H,Pmax_dBm,noise_pow);
SINR1=getSINR(H,W1,noise_pow)
[U_opt,p,p_ext,q,q_ext,n]=p1_iteration(H,Pmax_dBm,noise_pow);
W2=U_opt*diag(sqrt(p));
SINR_opt=getSINR(H,W2,noise_pow);

% data_former(H)
%[H_cell,q_cell,combine_cell]=data_former(N,K,100,d1_km,d2_km,Pmax_dBm,noise_pow)
[H_cell,q_cell,combine_cell]=data_former(N,K,25000,d1_km,d2_km,Pmax_dBm,noise_pow)

% [p_star,W_star]=p2_powmin(H,U_opt,SINR_opt,noise_pow);
% SINR_change=getSINR(H,W_star,noise_pow);
% [weight]=WMMSE(H,W,noise_pow)