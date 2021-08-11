function [U,p,p_ext,q,q_ext,n]=p1_iteration(H,Pmax_dBm,noise_pow)
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
        [ul_cpmat]=ul_cpmatrix(R,U,Pmax_dBm,noise_pow);
        [V_q,D_q]=eig(ul_cpmat);
        [lambda_maxval lambda_maxind]=max(real(diag(D_q)));
        
        q_ext=real(V_q(:,lambda_maxind)); %step 8
        q=q_ext(1:K)*1/q_ext(end); %?????????????????
        lambda_cell{n}=lambda_maxval;
        C{n}=1/lambda_cell{n}; %step 9
        if (n-1>0)&&((lambda_cell{n-1}-lambda_cell{n})<1e-5) %step 10
            break;
        end
    end
    %lambda_maxval
    [dl_cpmat]=dl_cpmatrix(R,U,Pmax_dBm,noise_pow);
    [V_p,D_p]=eig(dl_cpmat);
    [p_maxval p_maxind]=max(real(diag(D_p)));
    
    p_ext=V_p(:,p_maxind); %step 8
    p=p_ext(1:K)*1/p_ext(end); %?????????????????
    
end