function [p_star,W_star]=p2_powmin(H,W_opt,SINR_constraint,noise_pow) %W_opt=U_opt
    [N K]=size(H);
    l=ones(K,1);
    psi=zeros(K,K);
    for k=1:K
        for k1=1:K
            if k1==k
               psi(k,k1)=(1/SINR_constraint(k))*abs(H(:,k)'*W_opt(:,k))^2;
            else
               psi(k,k1)=-abs(H(:,k)'*W_opt(:,k1))^2;
            end
        end  
    end
    p_star=noise_pow*inv(psi)*l;
    W_star=W_opt*sqrt(diag(p_star));