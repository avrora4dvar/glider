function [rx] = beckman_haidvogel_haney(z)
%beckman_haidvolge_haney calculates beckman_haidvogel or Haney number. 
%   Calculates Beckman-Haidvogel (size(z,3)>1) or Haney number
%   (size(z,3)==1)
%   max_{i'j'kk'}=|(z_{ijk}-z_{i'j'k})-(z_{ijk-1}-z_{i'j'k-1)}|/...
%   |(h_{ijk}+h_{i'j'k})-(h_{ijk-1}+h_{i'j'k-1})| with maximum taken over adjacent
%   cells and all layers. The number must be small to avoid artifical horizontal
%   gradient terms caused by use sigma-grid. Non-active cells should be set
%   to NaN; 

if size(z,3)==1
    z=repmat(z,[1,1,2]); 
    z(:,:,1)=0; 
end

rx=zeros(size(z,1),size(z,2)); 
for k=2:size(z,3)
    rx(1:end-1,:)=nanmax(rx(1:end-1,:),...
        abs( (z(2:end,:,k)-z(1:end-1,:,k))+(z(2:end,:,k-1)-z(1:end-1,:,k-1)) )./...
        ( (z(2:end,:,k)+z(1:end-1,:,k))-(z(2:end,:,k-1)+z(1:end-1,:,k-1)) ) ); 
    rx(2:end,:)=nanmax(rx(2:end,:),...
        abs( (z(1:end-1,:,k)-z(2:end,:,k))+(z(1:end-1,:,k-1)-z(2:end,:,k-1)) )./...
        ( (z(1:end-1,:,k)+z(2:end,:,k))-(z(1:end-1,:,k-1)+z(2:end,:,k-1)) ) ); 
    rx(:,1:end-1)=nanmax(rx(:,1:end-1),...
        abs( (z(:,2:end,k)-z(:,1:end-1,k))+(z(:,2:end,k-1)-z(:,1:end-1,k-1)) )./...
        ( (z(:,2:end,k)+z(:,1:end-1,k))-(z(:,2:end,k-1)+z(:,1:end-1,k-1)) ) ); 
    rx(:,2:end)=nanmax(rx(:,2:end),...
        abs( (z(:,1:end-1,k)-z(:,2:end,k))+(z(:,1:end-1,k-1)-z(:,2:end,k-1)) )./...
        ( (z(:,1:end-1,k)+z(:,2:end,k))-(z(:,1:end-1,k-1)+z(:,2:end,k-1)) ) );

end

end

