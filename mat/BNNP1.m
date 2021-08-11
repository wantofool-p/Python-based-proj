function []=BNNP1(H)

    
    R_h=reshape(real(H),[],1);
    J_h=reshape(imag(H),[],1);
    
    input=[R_h,J_h]';
    
    height=N*K;
    width=2;
%     conv_layer=convolution2dLayer([height width],)
%      setwb
%      reluLayer
%      sigmoid
%        fullyConnectedLayer