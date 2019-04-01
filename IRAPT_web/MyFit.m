function [Res]=MyFit(X1,Y1,X2)
% Fast linear interpolation function
n1=length(X1);
n2=length(X2);
Res=zeros(1,n2);
Ind1=1;
for Ind2=1:n2
    while(Ind1<=n1 && X2(Ind2)>X1(Ind1))
        Ind1=Ind1+1;
    end
    if(Ind1>n1)
        break;
    end
    if(Ind1>1)
        Res(Ind2)=Y1(Ind1-1)+(Y1(Ind1)-Y1(Ind1-1))*(X2(Ind2)-X1(Ind1-1))/(X1(Ind1)-X1(Ind1-1));
    else
        Res(Ind2)=Y1(1);
    end
end
if(Ind1>n1)
    for Ind2=Ind2:n2
        Res(Ind2)=Y1(n1);
    end
end

end