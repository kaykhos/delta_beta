clear  %改变db和z，但是保证传输总长度不变
%H=xlsread('H-5site.xlsx');
%H=xlsread('H-qpr-original.xlsx');
H=xlsread('H-8site.xlsx');
H0=H(1:60,1:60);
H = H0;
%z=1;  % change delta beta at every 1mm
n=size(H,1);
%dbmax=0.0005*64;
%z=0.000001/dbmax/dbmax; 
dbmax=0.05;
z=1.5;

num=1;  %重复次数
meannum=30;%重复模拟100次取方差

%iternum=floor(16/z)+1;  %改变次数 总长3cm
iternum=16;
zl=mod(30,z);
Rho=zeros(n,iternum);  %终态  概率分布
rho_av=zeros(n,iternum);  %求和平均


R=zeros(iternum,7);
Ht=zeros(n,n);

chongfu = zeros(meannum,1);
for t=1:meannum
    gl = 0;
    xl = zeros(num,iternum); %记录中间态效率
    for calc=1:num
        R=dbmax.*rand([iternum,7]); %生成随机数
        %for i=1:iternum
        %    Ra=dbmax.*rand([1,n]);
        %    for j=i:n
        %        R(i,j)=Ra(1,j);
        %    end;
        %end;
        B=zeros(n,n);
        C=zeros(n,n);
        A=zeros(n,n);
        
        Psi=zeros(n,1);  %中间态
        % t=1/n;
        % t=sqrt(t);
        % for ii=1:n
        %    Psi(ii)=t;  %初态
        % end;
        Psi(6)=1;%从第六根注入
        %for i=1:1:n
        %   Psi(i)=1/sqrt(n);
        %end
        
        
        for ii=1:1:iternum
            for ge=1:1:7    %只有前7根波导考虑delta beta
                H(ge,ge)=H(ge,ge)-R(ii,ge);  %抽取R中的一行
            end
            for gep=1:1:6
                for geq=gep+1:1:7
                    if H(gep, geq)>0.15
%                         [gep,geq]
                        H(gep,geq)=sqrt(H(gep,geq)^2+(R(ii,gep)-R(ii,geq))^2/4);
                        H(geq,gep)=H(gep,geq);
                    end
                end
            end
            Ht=H;
            if ii==iternum
                B=H*(-1)*zl;
            else
                B=H*(-1)*z;
            end
            H=H0;
            C=B*sqrt(-1);
            A=expm(C);  %e^iHz
            Psi=A*Psi;
            Psi0=abs(Psi);
            for ge=1:1:n
                Rho(ge,ii)=Psi0(ge)*Psi0(ge);  %概率
                rho_av(ge,ii)=rho_av(ge,ii)+Rho(ge,ii);
            end
            for ge = 8:n
                gl=Psi0(ge)*Psi0(ge);  %概率
                xl(calc,ii)=xl(calc,ii)+gl;
            end
            
        end
        
    end
    rho_av=rho_av/num;
    
    out=zeros(iternum,1);
    pingjun=zeros(iternum,1);
    
    for x=iternum
        A = xl(:,x);
        pingjun(x)=mean(A);
        %plot(x,out);
        %hold on
    end  % line of process
    chongfu(t)=pingjun(iternum);
end
mean(chongfu)
std(chongfu)
