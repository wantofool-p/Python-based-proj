function [SINR]=p1_recovery(H,q)
    [N K]=size(H);
    n=0;
    q=zeros(K,1); %step 1
    
    R=cell(K,1);
    for i=1:K
        R{i}=H(:,i)*H(:,i)';
        R{i}=R{i}/noise_pow;    
    end         %step 2
    
    noise_pow=1; %step 3
    
    
    q_ext=[q;1];
    U=zeros(N,K);
    lambda_cell={};
    C={}
    while true %step 4
        n=n+1; %step 5
        for i=1:K      %step 6
            Q=eye(N);
            for k=1:K
                if k~=i
                Q=Q+q_ext(k)*R{k}; 
                end    
            end
            %target=R{i}*inv(Q)
            [Vtar,Dtar]=eig(R{i},Q);
            [maxval maxindex]=max(real(diag(Dtar)));
            Vtar_max=Vtar(:,maxindex);
            U(:,i)=Vtar_max; 
            U(:,i)=Vtar_max/norm(Vtar_max); %step 7
        end