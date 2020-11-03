
Hfilename='H-8site.xlsx';
Lfilename='L-8site.xlsx';

wgnum=27;  %7 sites plus 20 sink waveguides
H= zeros(wgnum,wgnum);
L= zeros(wgnum,wgnum);

Beta=0;
dbmax=1;  % vary dbmax from 0 to 1

for i =1:1:wgnum
 H(i,i)=Beta;
end

a= 0.0014;
C12=960*a;
C23=330*a;
C34=511*a; 
C45=766*a;
C46=142*a; 
C47=670*a;
C56=783*a;
C57=100*a;
C67=383*a;

Css=1; %coupling coefficient between sink waveguides


% Hamiltonian
H(1,2)= C12;
H(2,1)= C12; H(2,3)= C23;
H(3,2)= C23; H(3,4)= C34;  H(3,8)=Css;

H(4,3)=C34; H(4,5)= C45; H(4,7)= C47; H(4,6)= C46; 
H(5,4)=C45; H(5,6)= C56; H(5,7)= C57; 
H(6,5)=C56; H(6,7)= C67; H(6,4)= C46 ;
H(7,6)=C67; H(7,4)=C47; H(7,5)=C57; 
H(8,3)=Css; 

for i = 9:1:wgnum
    H(i-1,i)= Css ; 
    H(i,i-1)= Css; 
end 

% Lindblad
%dephasing part:
for i=1:7
    L(i,i)=sqrt(dbmax/2);  %delta beta follows a uniform distribution (0, dbmax), take the average deltabeta=dbmax/2 -> L=sqrt(delta beta)
end 

%decoherence part:
%Ceff^2=C0^2+(deltabeta/2)^2 -> (Ceff-C0)(Ceff+C0)=deltabeta^2/4, let delta C= Ceff-C0, Ceff+C0=2C0, ->delta C= deltabeta^2/8C0

% Delta beta for waveguide i and j both follow uniform distribution (0, dbmax), absolute value of (db1-db2) =dbmax/3 
for i=1:6
    for j=i+1:7
        if H(i,j)~=0
        L(i,j)=(dbmax/3)/sqrt(8)/sqrt(H(i,j));  % L=sqrt(delta C)=sqrt((db1-db2)^2/8C0)
        L(j,i)=L(i,j);
        end
    end 
end 

%L(8,3)=sqrt(Css); 
%L(3,8)= sqrt(Css);
L(8,3)=(dbmax/2)/sqrt(8)/sqrt(H(8,3));  % Waveguide 3 has delta beta and Waveguide 8 does not. 
L(3,8)=L(8,3);

Hred = H(1:8,1:8);
Lred = L(1:8,1:8);

%xlswrite('SevenSiteTransportEfficiency.xlsx', V1);
%xlswrite(Hfilename,H);
%xlswrite(Lfilename,L);