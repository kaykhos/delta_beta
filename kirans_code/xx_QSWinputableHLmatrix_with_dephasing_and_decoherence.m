Generated8SiteHL_sqrtC_decoherence_and_dephasing;
H = H(1:10, 1:10);
n=size(H,1);

timax=16;  % evolution time   evolution length 200mm
A=zeros(n,n);  % adjacency matrix
G=zeros(n,n);  % generator(google) matrix
I=eye(n,n);  % identity matrix
rho=zeros(n,n);  % density matrix
rho(6,6)=1;  % define the input node �ӵ�6��ע��
rhovd=reshape(rho,[n*n,1]); % vectorized rho
lk=zeros(n,n);  % Lk in the lindblad equation
lkcj=zeros(n,n);  % complex conjugate of lk
lindsum_dc=zeros(n.*n,n.*n);  % the sum used in lindblad equation for decoherence term
lindsum_dp=zeros(n.*n,n.*n);  % the sum used in lindblad equation for dephasing term

outcome=zeros(n,timax); % initialize outcome
for t=1:n
    for p=1:n
        if p~= t
        lk(t,p)=L(t,p);
        lkcj(t,p)=conj(L(t,p));
        lindsum_dc=lindsum_dc+kron(lkcj,lk)-1./2.*(kron(I,(lk'*lk))+kron((lk.'*lkcj),I));  % Page 14 equation 15
        lk(t,p)=0;
        lkcj(t,p)=0;
        end 
    end
end

for t=1:n
        lk(t,t)=L(t,t);
        lkcj(t,t)=conj(L(t,t));
        lindsum_dp=lindsum_dp+2.*(kron(lkcj,lk)-1./2.*(kron(I,(lk'*lk))+kron((lk.'*lkcj),I)));  % Page 14 equation 15
        lk(t,t)=0;
        lkcj(t,t)=0;
end

%lambda=-(1-omega).*sqrt(-1).*(kron(I,H)-kron(H.',I))+omega.*lindsum;
lambda=-sqrt(-1).*(kron(I,H)-kron(H.',I))+lindsum_dp+lindsum_dc;
for p=1:timax
    lambda0=lambda.*(p*10);
    rhovd1=expm(lambda0)*rhovd;
    for t=1:n
        outcome(t,p)=rhovd1(t+(t-1).*n,1);
    end
    outcome=abs(outcome);
end
out=zeros(timax,1);
for x=1:1:timax
       out(x)=1-(outcome(1,x)+outcome(2,x)+outcome(3,x)+outcome(4,x)+outcome(5,x)+outcome(6,x)+outcome(7,x)) %the probability at the sinks
        %plot(x,out);
        %hold on
end % line of process

figure(2);clf;shg
plot(out);
title('me')
xlabel('t')
outcome


