%% read the data

E_S=xlsread('dataS4.csv');%read energy levels of S
E_P=xlsread('dataP4.csv');%read energy levels of P
V_S=xlsread('dataSVe4.csv');%read eigenvectors of S
V_P=xlsread('dataPVe4.csv');%read eigenvectors of S

for k=1:3001 %build the matrix of energy
    for l=1:12     %lth eigenvector of S
        Norm=0;
        for j=1:12   %jth number
            Norm=Norm+V_S(k,(l-1)*12+j)^2;
        end
        for j=1:12  
            Vec_S(l,j,k) = V_S(k,(l-1)*12+j)/sqrt(Norm);
        end
    end
end

for k=1:3001 
    for l=1:24 %lth eigenvector of P
        Norm=0;
        for j=1:24
            Norm=Norm+V_P(k,(l-1)*24+j)^2;
        end
        for j=1:24
            Vec_P(l,j,k) = V_P(k,(l-1)*24+j)/sqrt(Norm);
        end
    end
end
%% Data initialization
B = 161;%magnetic field

%Unit conversion S:1-12  P:13-36
for l=1:12
    detuning(l)=1.012*10^(9)*E_S(B*10+1,l);
end
for l=13:36
    detuning(l)=25*10^(6)*E_P(B*10+1,l-12);
end
%detunings of imaging light
Detu = zeros(36,36);
%detunings of repumping light
Detu2 = zeros(36,36);

for l=1:12%imaging light
    for j=13:36
        Detu(l,j)=abs(detuning(l)-detuning(12)-detuning(j)+detuning(36));
        Detu(j,l)=Detu(l,j);
     end
end
for l=1:12
    for j=13:36
        Detu2(l,j)=abs(detuning(l)-detuning(10)-detuning(j)+detuning(34));
        Detu2(j,l)=Detu2(l,j);
    end
end

%transition coefficients
CG_o = zeros(36,36);%sigma+
CG_p = zeros(36,36);%pi
CG_n = zeros(36,36);%sigma-

s1=[1,2,4,6,8,10];%list
s2=[3,5,7,9,11,12];
p1=[13,14,16,19,23,27];
p2=[15,17,20,24,28,31];
p3=[18,21,25,29,32,34];
p4=[22,26,30,33,35,36];  
 %transition coefficients without IJ coupling 
for k=1:6
    CG_o(s1(k),p1(k))=1;%sigma+
    CG_o(p1(k),s1(k))=1;        
    CG_o(s2(k),p2(k))=sqrt(1/3);
    CG_o(p2(k),s2(k))=sqrt(1/3);  

    CG_p(s1(k),p2(k))=sqrt(2/3);    %pi
    CG_p(p2(k),s1(k))=sqrt(2/3);        
    CG_p(s2(k),p3(k))=sqrt(2/3);
    CG_p(p3(k),s2(k))=sqrt(2/3);  

    CG_n(s1(k),p3(k))=sqrt(1/3);%sigma-
    CG_n(p3(k),s1(k))=sqrt(1/3);        
    CG_n(s2(k),p4(k))=1;
    CG_n(p4(k),s2(k))=1;  
end
    %the eigenvector under a magnetic field
X=zeros(36,36);

for l=1:12
    for j=1:12
        X(l,j)=Vec_S(l,j,10*B+1);
    end
end
for l=13:36
    for j=13:36
        X(l,j)=Vec_P(l-12,j-12,10*B+1);
    end
end

%%

I=1;%intensity of imaging light
I2=0.4;%repumping light
Gamma_0=2*pi*6.0666*10^(6);% spontoneous scatering rate
Gamma=Gamma_0*I/2;%simulated radiation of imaging light
Gamma2=Gamma_0*I2/2;%simulated radiation of repumping light
De=X*(CG_o+CG_n+CG_p)*(X^-1);%transition matrix of spontoneous scatering
As=X*(CG_n)*(X^-1);          % simulated radiation of imaging light
As2=X*(0.707*CG_n+0.707*CG_o)*(X^-1);          %repumping
decay=zeros(36,36);%the rates including detunings
Asorb=zeros(36,36);
Asorb2=zeros(36,36);
    
for l=1:36
    for j=1:36%对应矩阵元
        decay(l,j)=Gamma_0*(De(l,j))^2;
        Asorb(l,j)=Gamma*(As(l,j))^2/(1+(4*pi*Detu(l,j)/Gamma_0)^2);      
        Asorb2(l,j)=Gamma2*(As2(l,j))^2/(1+(4*pi*Detu2(l,j)/Gamma_0)^2);%repump
    end
end
a=zeros(36,40001);%unmber distrubution
a(10,1)=1;
dt=0.000000001; %time intervals   1ns
N1(1)=0;%photon scattering numbers
for j=1:40000 %40 us
    N1(j+1)=N1(j); 
    for l=1:12 % S levels
        a(l,j+1)=a(l,j);
        for k=13:36
            a(l,j+1)=a(l,j+1)+decay(l,k)*a(k,j)*dt+Asorb(l,k)*a(k,j)*dt-Asorb(l,k)*a(l,j)*dt+Asorb2(l,k)*a(k,j)*dt-Asorb2(l,k)*a(l,j)*dt;
            N1(j+1)=N1(j+1)-Asorb(l,k)*a(k,j)*dt+Asorb(l,k)*a(l,j)*dt;
        end
    end
    for l=13:36 % P levels
        a(l,j+1)=a(l,j);
        for k=1:12
            a(l,j+1)=a(l,j+1)-decay(k,l)*a(l,j)*dt+Asorb(l,k)*a(k,j)*dt-Asorb(l,k)*a(l,j)*dt+Asorb2(l,k)*a(k,j)*dt-Asorb2(l,k)*a(l,j)*dt;
        end        
    end
end

%% pictures 
figure(4)
plot(N1);



