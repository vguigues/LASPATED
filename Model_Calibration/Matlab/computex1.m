

function [y]=computex1(i,j,type,bs,rs)

piapp=3.14159;
y=0;
if (type==1)
    %red
    for k=1:10
        int3=ceil((k/5)*(i-1));
        int4=floor((k/5)*i);
        r1=(5/k)*int3;
        if (r1<=i)
            pw1=(-1)^(int3);
            int1=(5/(piapp*k))*pw1*(cos(piapp*int3)-cos(piapp*k*((i-1)/5)));
            int1=int1+(10/(piapp*k))*(int4-int3)-(5/(piapp*k))*((-1)^(int4))*(cos(piapp*k*(i/5))-cos(piapp*int4));
        else
            int1=(5/piapp*k)*((-1)^(int4))*(cos(piapp*k*((i-1)/5))-cos(piapp*k*i/5));    
        end
        int3=ceil((k/5)*(j-1));
        int4=floor((k/5)*j);
        r1=(5/k)*int3;
        if (r1<=j)
            pw1=(-1)^(int3);
            int2=(5/(piapp*k))*pw1*(cos(piapp*int3)-cos(piapp*k*((j-1)/5)));
            int2=int2+(10/(piapp*k))*(int4-int3)-(5/(piapp*k))*((-1)^(int4))*(cos(piapp*k*(j/5))-cos(piapp*int4));
        else
            int2=(5/piapp*k)*((-1)^(int4))*(cos(piapp*k*((j-1)/5))-cos(piapp*k*j/5));    
        end
        y=y+(1/(2^k))*((rs(2*k-1)+1)*int1+(rs(2*k)+1)*int2);
    end
else
    %blue
    for k=1:10
        int3=ceil((k/5)*(i-1));
        int4=floor((k/5)*i);
        r1=(5/k)*int3;
        if (r1<=i)
            pw1=(-1)^(int3);
            int1=(5/(piapp*k))*pw1*(cos(piapp*int3)-cos(piapp*k*((i-1)/5)));
            int1=int1+(10/(piapp*k))*(int4-int3)-(5/(piapp*k))*((-1)^(int4))*(cos(piapp*k*(i/5))-cos(piapp*int4));
        else
            int1=(5/piapp*k)*((-1)^(int4))*(cos(piapp*k*((i-1)/5))-cos(piapp*k*i/5));    
        end
        int3=ceil((k/5)*(j-1));
        int4=floor((k/5)*j);
        r1=(5/k)*int3;
        if (r1<=j)
            pw1=(-1)^(int3);
            int2=(5/(piapp*k))*pw1*(cos(piapp*int3)-cos(piapp*k*((j-1)/5)));
            int2=int2+(10/(piapp*k))*(int4-int3)-(5/(piapp*k))*((-1)^(int4))*(cos(piapp*k*(j/5))-cos(piapp*int4));
        else
            int2=(5/piapp*k)*((-1)^(int4))*(cos(piapp*k*((j-1)/5))-cos(piapp*k*j/5));    
        end
        y=y+(1/(2^k))*((bs(2*k-1))*int1+(bs(2*k))*int2);
    end    
end
      