function [H,distance,small_scaling]=hk_producer(d1_km,d2_km,N,K)
    distance=(d2_km-d1_km).*rand(K,1)+d1_km; %UE position in disc with two boundries
    path_loss=128.1+37.6*log10(distance);%loss abs
    H=zeros(N,K);
    for user=1:K
        small_scaling=sqrt(1/2)*(randn(N,1)+randn(N,1)*1i);
        H(:,user)=sqrt(10.^(-path_loss(user)/10)).*small_scaling; %channel feature
    end
end