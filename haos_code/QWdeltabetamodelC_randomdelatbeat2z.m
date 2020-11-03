clear  %�ı�db��z�����Ǳ�֤�����ܳ��Ȳ���
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

num=1;  %�ظ�����
meannum=30;%�ظ�ģ��100��ȡ����

%iternum=floor(16/z)+1;  %�ı���� �ܳ�3cm
iternum=16;
zl=mod(30,z);
Rho=zeros(n,iternum);  %��̬  ���ʷֲ�
rho_av=zeros(n,iternum);  %���ƽ��


R=zeros(iternum,7);
Ht=zeros(n,n);

chongfu = zeros(meannum,1);
for t=1:meannum
    gl = 0;
    xl = zeros(num,iternum); %��¼�м�̬Ч��
    for calc=1:num
        R=dbmax.*rand([iternum,7]); %���������
        %for i=1:iternum
        %    Ra=dbmax.*rand([1,n]);
        %    for j=i:n
        %        R(i,j)=Ra(1,j);
        %    end;
        %end;
        B=zeros(n,n);
        C=zeros(n,n);
        A=zeros(n,n);
        
        Psi=zeros(n,1);  %�м�̬
        % t=1/n;
        % t=sqrt(t);
        % for ii=1:n
        %    Psi(ii)=t;  %��̬
        % end;
        Psi(6)=1;%�ӵ�����ע��
        %for i=1:1:n
        %   Psi(i)=1/sqrt(n);
        %end
        
        
        for ii=1:1:iternum
            for ge=1:1:7    %ֻ��ǰ7����������delta beta
                H(ge,ge)=H(ge,ge)-R(ii,ge);  %��ȡR�е�һ��
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
                Rho(ge,ii)=Psi0(ge)*Psi0(ge);  %����
                rho_av(ge,ii)=rho_av(ge,ii)+Rho(ge,ii);
            end
            for ge = 8:n
                gl=Psi0(ge)*Psi0(ge);  %����
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
