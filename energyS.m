energy=zeros(12,3001);
vector=zeros(12,12,3001);
HC=zeros(12,12);%coupling
H=zeros(12,12);
HB=zeros(12,12);%magnetic field
%submatrix
S2=zeros(2,2);
S3=zeros(2,2);
S4=zeros(2,2);
S5=zeros(2,2);
S6=zeros(2,2);

%Hamitonian of IJ coupling
HC(1,1)=5/4;%m_F=+-3
HC(12,12)=5/4;

HC(2,2)=3/4;%m_F=2
HC(3,3)=-5/4;
HC(2,3)=sqrt(5)/2;
HC(3,2)=sqrt(5)/2;

HC(4,4)=1/4;%m_F=1
HC(5,5)=-3/4;
HC(4,5)=sqrt(2);
HC(5,4)=sqrt(2);

HC(6,6)=-1/4;%m_F=0
HC(7,7)=-1/4;
HC(6,7)=3/2;
HC(7,6)=3/2;

HC(8,8)=-3/4;%m_F=-1
HC(9,9)=1/4;
HC(8,9)=sqrt(2);
HC(9,8)=sqrt(2);

HC(10,10)=-5/4;%m_F=-2
HC(11,11)=3/4;
HC(10,11)=sqrt(5)/2;
HC(11,10)=sqrt(5)/2;

for i=1:3001
    B = 0.1*i-0.1;%magnetic field
    c1 = 0.0027691*B;
    c2 = 0.0000004061*B;%Hamitonian of magnetic field
    HB(1,1) = 0.5*c1+2.5*c2;
    HB(2,2) = 0.5*c1+1.5*c2;  
    HB(3,3) = -0.5*c1+2.5*c2;  
    HB(4,4) = 0.5*c1+0.5*c2;   
    HB(5,5) = -0.5*c1+1.5*c2;   
    HB(6,6) = 0.5*c1-0.5*c2;    
    HB(7,7) = -0.5*c1+0.5*c2;   
    HB(8,8) = 0.5*c1-1.5*c2;   
    HB(9,9) = -0.5*c1-0.5*c2;   
    HB(10,10) = 0.5*c1-2.5*c2;   
    HB(11,11) = -0.5*c1-1.5*c2;   
    HB(12,12) = -0.5*c1-2.5*c2;  % Hamitonian of magnetic field
    
    H=HC+HB;
    
    energy(1,i)=H(1,1);%mF=+3
    vector(1,1,i)=1;
    energy(12,i)=H(12,12);%mF=-3
    vector(12,12,i)=1;
 
    for j=1:2%mF=+-2
        for k=1:2
            S2(j,k)=H(j+1,k+1);
            S3(j,k)=H(j+3,k+3);
            S4(j,k)=H(j+5,k+5);
            S5(j,k)=H(j+7,k+7);
            S6(j,k)=H(j+9,k+9);
        end        
    end
    [V2,D2] = eig(S2);
    [d2,ind2] = sort(diag(D2));
    [V3,D3] = eig(S3);
    [d3,ind3] = sort(diag(D3));
    [V4,D4] = eig(S4);
    [d4,ind4] = sort(diag(D4));
    [V5,D5] = eig(S5);
    [d5,ind5] = sort(diag(D5));
    [V6,D6] = eig(S6);
    [d6,ind6] = sort(diag(D6));

    for j=1:2
        energy(j+1,i)=d2(j);
        energy(j+3,i)=d3(j);
        energy(j+5,i)=d4(j);
        energy(j+7,i)=d5(j);
        energy(j+9,i)=d6(j);%sort the eigenvalues
        for k=1:2
            vector(j+1,k+1,i)=V2(j,ind2(k));
            vector(j+3,k+3,i)=V3(j,ind3(k));
            vector(j+5,k+5,i)=V4(j,ind4(k));
            vector(j+7,k+7,i)=V5(j,ind5(k));
            vector(j+9,k+9,i)=V6(j,ind6(k));%sort the eigenvectors
        end
    end
end  
figure(1)
plot(energy');

V=zeros(3001,144);
for i=1:3001
    for j=1:12
        for k=1:12
            V(i,(j-1)*12+k)=vector(k,j,i);
        end
    end
end
%% export the data    
writematrix(energy','dataS4.csv'); 
writematrix(V,'dataSVe4.csv');
    