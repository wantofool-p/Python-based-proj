function [C]=dl_cpmatrix(R,W,Pmax_dBm,noise_pow)
    [N K]=size(W);
    Pmax_linear=10^((Pmax_dBm-30)/10);
    
    l=ones(K,1);
    v=noise_pow*l; 
    U=zeros(K,K);
    D=zeros(K,K);
    for k=1:K
        for k1=1:K
            if k1~=k
               U(k,k1)=W(:,k1)'*R{k}*W(:,k1);
               D(k,k1)=0;
            else
               U(k,k1)=0;
               D(k,k1)=1/(W(:,k1)'*R{k}*W(:,k1));
            end
        end  
    end
    C=[D*U D*v; (1/Pmax_linear)*l'*D*U (1/Pmax_linear)*l'*D*v];
end