energy=zeros(24,3001);
vector=zeros(24,24,3001);
HC=zeros(24,24);%coupling
H=zeros(24,24);
HB=zeros(24,24);%magnetic field
%submatrix
S2=zeros(2,2);
S3=zeros(3,3);
S4=zeros(4,4);
S5=zeros(4,4);
S6=zeros(4,4);
S7=zeros(3,3);
S8=zeros(2,2);

%Hamitonian of IJ coupling
HC(1,1)=15/4;%m_F=+-4
HC(24,24)=15/4;

HC(2,2)=9/4;%m_F=3
HC(3,3)=5/4;
HC(2,3)=sqrt(15)/2;
HC(3,2)=sqrt(15)/2;

HC(23,23)=9/4;%m_F=-3
HC(22,22)=5/4;
HC(23,22)=sqrt(15)/2;
HC(22,23)=sqrt(15)/2;

HC(4,4)=3/4;%m_F=2
HC(5,5)=3/4;
HC(6,6)=-5/4;
HC(4,5)=sqrt(6);
HC(5,4)=sqrt(6);
HC(6,5)=sqrt(5);
HC(5,6)=sqrt(5);

HC(21,21)=3/4;%m_F=-2
HC(20,20)=3/4;
HC(19,19)=-5/4;
HC(21,20)=sqrt(6);
HC(20,21)=sqrt(6);
HC(19,20)=sqrt(5);
HC(20,19)=sqrt(5);

HC(7,7)=-3/4;%m_F=1
HC(8,8)=1/4;
HC(9,9)=-3/4;
HC(10,10)=-15/4;
HC(7,8)=sqrt(27)/2;
HC(8,7)=sqrt(27)/2;
HC(9,8)=sqrt(8);
HC(8,9)=sqrt(8);
HC(9,10)=sqrt(15)/2;
HC(10,9)=sqrt(15)/2;

HC(18,18)=-3/4;%m_F=-1
HC(17,17)=1/4;
HC(16,16)=-3/4;
HC(15,15)=-15/4;
HC(17,18)=sqrt(27)/2;
HC(18,17)=sqrt(27)/2;
HC(16,17)=sqrt(8);
HC(17,16)=sqrt(8);
HC(16,15)=sqrt(15)/2;
HC(15,16)=sqrt(15)/2;

HC(11,11)=-9/4;%m_F=0;
HC(12,12)=-1/4;
HC(13,13)=-1/4;
HC(14,14)=-9/4;
HC(11,12)=sqrt(6);
HC(12,11)=sqrt(6);
HC(12,13)=3;
HC(13,12)=3;
HC(14,13)=sqrt(6);
HC(13,14)=sqrt(6);

for i=1:3001
    B = 0.1*i-0.1;
    c1 = 0.0748067116403314*B;%Hamitonian of magnetic field
    c2 = -0.000021966242806066916*B;
    HB(1,1) = 1.5*c1+2.5*c2;
    HB(2,2) = 1.5*c1+1.5*c2;  
    HB(3,3) = 0.5*c1+2.5*c2;  
    HB(4,4) = 1.5*c1+0.5*c2;   
    HB(5,5) = 0.5*c1+1.5*c2;   
    HB(6,6) = -0.5*c1+2.5*c2;    
    HB(7,7) = 1.5*c1-0.5*c2;   
    HB(8,8) = 0.5*c1+0.5*c2;   
    HB(9,9) = -0.5*c1+1.5*c2;   
    HB(10,10) = -1.5*c1+2.5*c2;   
    HB(11,11) = 1.5*c1-1.5*c2;   
    HB(12,12) = 0.5*c1-0.5*c2;   
    HB(13,13) = -0.5*c1+0.5*c2;   
    HB(14,14) = -1.5*c1+1.5*c2;
    HB(24,24) = -1.5*c1-2.5*c2;
    HB(23,23) = -1.5*c1-1.5*c2;  
    HB(22,22) = -0.5*c1-2.5*c2;  
    HB(21,21) = -1.5*c1-0.5*c2;   
    HB(20,20) = -0.5*c1-1.5*c2;   
    HB(19,19) = 0.5*c1-2.5*c2;    
    HB(18,18) = -1.5*c1+0.5*c2;   
    HB(17,17) = -0.5*c1-0.5*c2;   
    HB(16,16) = 0.5*c1-1.5*c2;   
    HB(15,15) = 1.5*c1-2.5*c2;
    
    H=1.025*HC+0.05*HC*HC+HB;%Hamitonian including electric quadrupole, magnetic field, IJ coupling 
    
    energy(1,i)=H(1,1);%mF=+4
    vector(1,1,i)=1;
    energy(24,i)=H(24,24);%mF=-4
    vector(24,24,i)=1;
 %construct submatrix 
    for j=1:2%mF=+-3
        for k=1:2
            S2(j,k)=H(j+1,k+1);
            S8(j,k)=H(21+j,21+k);
        end        
    end
    [V2,D2] = eig(S2);
    [d2,ind2] = sort(diag(D2));
    [V8,D8] = eig(S8);
    [d8,ind8] = sort(diag(D8));
    for j=1:2
        energy(j+1,i)=d2(j);
        energy(21+j,i)=d8(j);%sort the eigenvalues
        for k=1:2
            vector(j+1,k+1,i)=V2(j,ind2(k));
            vector(21+j,21+k,i)=V8(j,ind8(k));%sort the eigenvectors
        end
    end
    
    for j=1:3%mF=+-2
        for k=1:3
            S3(j,k)=H(j+3,k+3);
            S7(j,k)=H(18+j,18+k);
        end        
    end
    [V3,D3] = eig(S3);;
    [V7,D7] = eig(S7);
    [d3,ind3] = sort(diag(D3))
    [d7,ind7] = sort(diag(D7));
    for j=1:3
        energy(j+3,i)=d3(j);
        energy(18+j,i)=d7(j);
        for k=1:3
            vector(j+3,k+3,i)=V3(j,ind3(k));
            vector(18+j,18+k,i)=V7(j,ind7(k));
        end
    end    
    
    for j=1:4%mF=+-1 0
        for k=1:4
            S4(j,k)=H(j+6,k+6);
            S5(j,k)=H(j+10,k+10);
            S6(j,k)=H(14+j,14+k);
        end        
    end
    [V4,D4] = eig(S4);
    [d4,ind4] = sort(diag(D4));
    [V5,D5] = eig(S5);
    [d5,ind5] = sort(diag(D5));
    [V6,D6] = eig(S6);
    [d6,ind6] = sort(diag(D6));    
    for j=1:4
        energy(j+6,i)=d4(j);
        energy(j+10,i)=d5(j);
        energy(j+14,i)=d6(j);
        for k=1:4
            vector(j+6,k+6,i)=V4(j,ind4(k));
            vector(j+10,k+10,i)=V5(j,ind5(k));
            vector(14+j,14+k,i)=V6(j,ind6(k));
        end
    end
end
figure(1)
plot(energy');

V=zeros(3001,576);
for i=1:3001
    for j=1:24
        for k=1:24
            V(i,(j-1)*24+k)=vector(k,j,i);
        end
    end
end
%%    export the data
writematrix(energy','dataP4.csv'); 
writematrix(V,'dataPVe4.csv');
    
    
    

