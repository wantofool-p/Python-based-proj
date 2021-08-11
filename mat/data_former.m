function [H_cell,q_cell,combine_cell]=data_former(N,K,num_samples,d1_km,d2_km,Pmax_dBm,noise_pow)
    H_cell=cell(num_samples,1);
    q_cell=cell(num_samples,1);
    combine_cell=cell(num_samples,4);
    combine_cell(1,:)={'No','Real(h)','Imag(h)','q'} %header
    for t=1:num_samples
        [H,d,n_complex]=hk_producer(d1_km,d2_km,N,K);
        h_bar=zeros(K,N)         %6*4-->4*6
        for k=1:K
            h_bar(k,:)=transpose(H(:,k));
        end
        H_cell{t}=H;
        R_h=reshape(real(h_bar)/sqrt(noise_pow),[],1); %24*1
        J_h=reshape(imag(h_bar)/sqrt(noise_pow),[],1);
        [U_opt,p,p_ext,q,q_ext,n]=p1_iteration(H,Pmax_dBm,noise_pow);
        q_cell{t}=q;
        combine_cell{t+1,1}=t;
        combine_cell{t+1,2}=R_h;
        combine_cell{t+1,3}=J_h;
        combine_cell{t+1,4}=q;
    end    
    T=cell2table(combine_cell(2:end,:),'VariableNames',combine_cell(1,:));
    writetable(T,'MyData_25000.csv')
    %writetable(T,'MyData_200.csv')
end