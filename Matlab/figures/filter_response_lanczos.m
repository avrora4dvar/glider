function [H] = filter_response_lanczos(f,dt,fcut,order)
%FILTER_RESPONSE_LANCZOS calculates the filter response function for the
%cosine-Lanczos filter
% [H]=filter_response_lanczos(f,dt,fcut) calculates the filter response function in frequency domain at
% frequencies f for for Lanczos filter with time-step dt and cut-off
% frequency fcut

step=1e-3;
fIn=[0:step:.5-step]; 
t=[0:1/step-1]; 

H=nan(size(fIn)); 
for k=1:length(fIn); 
signalIn=cos(2*pi*t*fIn(k)); 
Fin=fft(signalIn)/length(signalIn);
Fin=Fin(1:length(fIn)); 
[dummy,imax]=max(Fin); 

signalOut=filter_lanczos(signalIn,1/dt,fcut,order,2); 
Fout=fft(signalOut)/length(signalOut);

H(imax)=Fout(imax)/Fin(imax); 
end


H=interp1(fIn,H,abs(f*dt)); 

end

