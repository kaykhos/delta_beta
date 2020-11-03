%%
%H=xlsread('H-5site.xlsx');
%H=xlsread('H-qpr-original.xlsx');
%H=xlsread('H-8site.xlsx');
H0=H; 

n=size(H,1);
dbmax=1;
num=200;  %Monte carlo repeat times

iternum=50;  % number of segments
z=1;  % change delta beta at every 2mm
Rho=zeros(n,iternum);  %probability distribution
sum=zeros(n,iternum);  %


R=zeros(iternum,7);
Ht=zeros(n,n);
for calc=1:num
    R=dbmax.*rand([iternum,7]); % random number generation
    B=zeros(n,n);
    C=zeros(n,n);
    A=zeros(n,n);

    Psi=zeros(n,1);  %
 
    Psi(6)=1;% input from waveguide no. 6
   

    for ii=1:1:iternum
        for ge=1:1:7    % only the 7 sites consider delta beta
            H(ge,ge)=H(ge,ge)-R(iternum,ge);  %
        end 
%         for gep=1:1:6
%             for geq=gep+1:1:7
%                 if H(gep, geq)>0.15
%                     H(gep,geq)=sqrt(H(gep,geq)^2+(R(iternum,gep)-R(iternum,geq))^2/4);
%                     H(geq,gep)=H(gep,geq);
%                 end
%             end 
%         end

            Ht=H;
            B=H*(-1)*z;
            H=H0; 
            C=B*sqrt(-1);
            A=expm(C);  %e^iHz
            Psi=A*Psi; 
            Psi0=abs(Psi);
     for ge=1:1:n
            Rho(ge,ii)=Psi0(ge)*Psi0(ge);  %
            sum(ge,ii)=sum(ge,ii)+Rho(ge,ii);
     end 
    end
   
end
sum=sum/num;
disp('done')
figure(10)
clf


pcolor(sum)
%%
if 0
    figure(1)
    subplot(2, 2 ,1)
    pcolor(Rho0)
    xlabel('segment #')
    ylabel('waveguide #')
    title('\Delta \beta = 0')
    axis square

    subplot(2, 2 ,2)
    pcolor(Rho1)
    xlabel('segment #')
    ylabel('waveguide #')
    title('\Delta \beta = 0.1')
    axis square

    subplot(2, 2, 3)
    pcolor(Rho2)
    xlabel('segment #')
    ylabel('waveguide #')
    title('\Delta \beta = 1')
    axis square

    subplot(2, 2, 4)
    pcolor(Rho3)
    xlabel('segment #')
    ylabel('waveguide #')
    title('\Delta \beta = 5')
    axis square
end

%%
% lets look at eiven vectors
[vecs, vals] = eigs(H + 5*eye(27), 27);
legendCell = cellstr(num2str(diag(vals), 'E=%-d'));
c = parula(27);

figure(15);clf;shg
subplot(2, 2, 1);hold on
plot(vecs(:,1:7))
plot([6,6], [-0.6, 0.6], 'r')
legend(legendCell(1:7))
title('Eigenstate of the (un)perturbed Hamiltonian')
xlabel('site #')
ylabel('amplitude')

subplot(2, 2, 2);hold on
plot(vecs(:,8:14))
plot([6,6], [-0.6, 0.6], 'r')
legend(legendCell(8:14))
title('Eigenstate of the (un)perturbed Hamiltonian')
xlabel('site #')
ylabel('amplitude')

subplot(2, 2, 3);hold on
plot(vecs(:,14:21))
plot([6,6], [-0.6, 0.6], 'r')
legend(legendCell(14:21))
title('Eigenstate of the (un)perturbed Hamiltonian')
xlabel('site #')
ylabel('amplitude')

subplot(2, 2, 4);hold on
plot(vecs(:,22:27))
plot([6,6], [-0.6, 0.6], 'r')
legend(legendCell(22:27))
title('Eigenstate of the (un)perturbed Hamiltonian')
xlabel('site #')
ylabel('amplitude')
%%

out=zeros(iternum,1);
for x=1:1:iternum
    %out(x)=sum(3,x);
    for ge=8:1:n
        out(x)=out(x)+sum(ge,x); % transfer efficiency to sinks
    end
   
	%plot(x,out);
	%hold on
end; % line of process
figure(1);clf;shg;
plot(out); 
title('monte carlo')
title('num')

Hred = H(1:7,1:7)
Lred = L(1:7,1:7)