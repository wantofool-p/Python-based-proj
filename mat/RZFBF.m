function [W,D]=RZFBF(H,Pmax_dBm,noise_pow)
     Pmax_linear=10^((Pmax_dBm-30)/10);
     [N K]=size(H);
     M=(pinv(H*(H)'+noise_pow.*eye(N))*H)
     W=sqrt(Pmax_linear)*M/norm(M,'fro');
     D=H'*W; %recover check
end