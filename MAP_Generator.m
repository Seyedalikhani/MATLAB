function a = MAP_Generator(SubCarrier_Mat)



% function y = truncated_sinc(t, fsamp, M)
M=64;
Num_samples=M;
Dist=randi([0 M*(M+1)/2],1,Num_samples);
MAP_Message=[];

MAP_Inerval=[];
for p=1:M
    MAP_Inerval=[MAP_Inerval; p (p-1)*(p)/2+1 p*(p+1)/2];
end



for g=1:Num_samples
    for j=1:M
        if (Dist(g)>=MAP_Inerval(j,2) && Dist(g)<=MAP_Inerval(j,3))
            MAP_Message=[MAP_Message MAP_Inerval(j,1)] ;
            break
        end
    end
end






a=SubCarrier_Mat(MAP_Message);


end
% 
% M1=sqrt(M);
% SubCarrier_Mat=zeros(M1,M1);
% Message_Mat=zeros(M1,M1);
% for m=1:M1
%     for n=1:M1
%         SubCarrier_Mat(n,m)=-M1+1+2*(m-1)-(-M1+1+2*(n-1))*1i;
%         Message_Mat(n,m)=n+8*(m-1);
%     end
% end
% 
% 
% 
% 
% a=SubCarrier_Mat(MAP_Message)+(-0.4+rand(1)*0.8+(-0.4+rand(1)*0.8)*1i);
% plot(real(a),imag(a),'b*')
% hold on
% grid on
% title('A noisy diagram of MAP Message')
% 
% 
