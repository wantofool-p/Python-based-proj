function [SINR]=getSINR(H,W,noise_pow)
    [N K]=size(H);
    SINR=zeros(K,1);
    for k=1:K
        signal=(abs((H(:,k))'*W(:,k)))^2;
        interference=0;
        for k1=1:K
            if k1~=k
            interference=interference+abs((H(:,k))'*W(:,k1))^2;
            end
        end
    SINR(k)=signal/(interference+noise_pow);
    end
end