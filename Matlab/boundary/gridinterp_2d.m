function [a]=gridinterp_2d(XX,YY,A,in0,iclosest,xb,yb)

A1=A;
if ~isempty(in0)
 A1(in0)=A1(iclosest); % since in0, iclosest are 1d, operate on a 2d array
end;
a=interp2(YY,XX,A1,yb,xb);
