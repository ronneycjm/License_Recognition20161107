function [n1,n2]=Find2Min(s)

[~,n1]=min(s);
s(n1)=s(n1)+10000;
[~,n2]=min(s);