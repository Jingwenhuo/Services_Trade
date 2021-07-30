%  Copyright ?  2018 Meng Jing All rights reserved

clear all;  %清变量
clc;    %清屏
tic;     %计时器，开始

Deflator=xlsread('C:\Users\14718\Desktop\IO\MRIO\priceindex_2010-2018.xlsx','Deflator');
Country_map=xlsread('C:\Users\14718\Desktop\IO\MRIO\priceindex_2010-2018.xlsx','Country_map');
load('C:\Users\14718\Desktop\IO\MRIO\CarbonData\MRIO_CO2_new.mat');  
%MRIO2018=xlsread('F:\Income groups_consumption categories\ADB2018\ADB_MRIO_2018.xlsx');
%%Energy/CO2
CO2_energy1=cell2mat(struct2cell(load('C:\Users\14718\Desktop\IO\MRIO\EnergyData\CO2_Energy_new2018.mat')));
CO2_energy1=CO2_energy1./1000;
original=CO2_energy1;
OCO2(:,1)=sum(original(:,1:3,1),2);
OCO2(:,2)=sum(original(:,1:3,2),2);
OCO2(:,3)=sum(original(:,1:3,3),2);
OCO2(:,4)=sum(original(:,1:3,4),2);
OCO2(:,5)=sum(original(:,1:3,5),2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
OCO2(:,6)=sum(original(:,1:3,6),2);
OCO2(:,7)=sum(original(:,1:3,7),2);
OCO2(:,8)=sum(original(:,1:3,8),2);
OCO2(:,9)=sum(original(:,1:3,9),2);

MRIO_energy=cell2mat(struct2cell(load('C:\Users\14718\Desktop\IO\MRIO\EnergyData\MRIO_Energy_2010-2018.mat')));
%MRIO_energy2017=cell2mat(struct2cell(load('C:\Users\HP\Desktop\IO\MRIO\EnergyData\MRIO_Energy_2010-2017.mat')));
%MRIO_energy=load('C:\Users\14718\Desktop\IO\MRIO\EnergyData\MRIO_Energy_2017.mat');
% 
energy_new2(:,:,1)=cell2mat(struct2cell(MRIO_energy(1)));
energy_new2(:,:,2)=cell2mat(struct2cell(MRIO_energy(2)));
energy_new2(:,:,3)=cell2mat(struct2cell(MRIO_energy(3)));
energy_new2(:,:,4)=cell2mat(struct2cell(MRIO_energy(4)));
energy_new2(:,:,5)=cell2mat(struct2cell(MRIO_energy(5)));
energy_new2(:,:,6)=cell2mat(struct2cell(MRIO_energy(6)));
energy_new2(:,:,7)=cell2mat(struct2cell(MRIO_energy(7)));
energy_new2(:,:,8)=cell2mat(struct2cell(MRIO_energy(8)));
energy_new2(:,:,9)=cell2mat(struct2cell(MRIO_energy(9)));

energy_new2=abs(energy_new2);%全部取正

emissionintensity=zeros(35*59,35*59*4,9);
share=zeros(4*35*59,35*59,9);
energyintensity=zeros(35*59,1,9);
%CO2_energy=permute(CO2_energy,[3,2,1]);% coal，oil, gas
%energy=permute(energy_new2,[4,3,2,1]); % coal，oil, gas， other

for k=1:9
    a=CO2_energy1(:,1:3,k);
    CO2_energy(:,:,k)= a';
    energy(:,:,k)= energy_new2(:,:,k);
end
energy_=zeros(4,35*59,9);

for k=1:9
    energy_(:,:,k)=energy(:,:,k)';
end  

abnormal=0;
for j=1:3
    for i=1:35*59
         for k=1:9
             if energy_(j,i,k)==0 && CO2_energy(j,i,k)~=0
                 abnormal=abnormal+1;
                 CO2_energy(j,i,k)=0;
             end
         end
    end
end

for k=1:9
    a=CO2_energy(1:3,:,k);
    CO2_energy1(:,1:3,k)= a';
    energy(:,:,k)= energy_new2(:,:,k);
end

%% column: 59*35+59*5+1; 2361
%% row: 59*35+1+6+1+1; 2084
CO2=zeros(35*59,9);%总排放就是三种能源的和,不要用之前那个总排放,不要waste
p1=CO2_energy1(:,1:4,1);
p2=sum(CO2_energy1(:,1:3,1),2);
CO2(:,1)=sum(CO2_energy1(:,1:3,1),2);
CO2(:,2)=sum(CO2_energy1(:,1:3,2),2);
CO2(:,3)=sum(CO2_energy1(:,1:3,3),2);
CO2(:,4)=sum(CO2_energy1(:,1:3,4),2);
CO2(:,5)=sum(CO2_energy1(:,1:3,5),2)
CO2(:,6)=sum(CO2_energy1(:,1:3,6),2);
CO2(:,7)=sum(CO2_energy1(:,1:3,7),2);
CO2(:,8)=sum(CO2_energy1(:,1:3,8),2);
CO2(:,9)=sum(CO2_energy1(:,1:3,9),2);
P=sum(CO2(:,:));

production_based=[sum(CO2(7*35+1:8*35,:))
    sum(CO2(10*35+1:11*35,:))
    sum(CO2(16*35+1:17*35,:))
    sum(CO2(21*35+1:22*35,:))
    sum(CO2(24*35+1:25*35,:))
    sum(CO2(42*35+1:43*35,:))
    sum(CO2(:,:))];
    
Energy_zero=(sum(OCO2)-sum(CO2))./sum(OCO2);%CO2异常值占排放的多少

M=zeros(59*35+1+6+1,59*35+59*5+1,9);
M(:,:,1)=CO2_2010.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,2)=CO2_2011.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,3)=CO2_2012.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,4)=CO2_2013.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,5)=CO2_2014.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,6)=CO2_2015.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,7)=CO2_2016.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,8)=CO2_2017.MRIO_CO2(1:59*35+1+6+1,:);
M(:,:,9)=CO2_2018.MRIO_CO2(1:59*35+1+6+1,:);

V(:,:,1)=sum(CO2_2010.MRIO_CO2(2067:2072,1:35*59));
V(:,:,2)=sum(CO2_2011.MRIO_CO2(2067:2072,1:35*59));
V(:,:,3)=sum(CO2_2012.MRIO_CO2(2067:2072,1:35*59));
V(:,:,4)=sum(CO2_2013.MRIO_CO2(2067:2072,1:35*59));
V(:,:,5)=sum(CO2_2014.MRIO_CO2(2067:2072,1:35*59));
V(:,:,6)=sum(CO2_2015.MRIO_CO2(2067:2072,1:35*59));
V(:,:,7)=sum(CO2_2016.MRIO_CO2(2067:2072,1:35*59));
V(:,:,8)=sum(CO2_2017.MRIO_CO2(2067:2072,1:35*59));
V(:,:,9)=sum(CO2_2018.MRIO_CO2(2067:2072,1:35*59));
for i=1:59
    for k=1:9
     V(1,(i-1)*35+1:i*35,k)=V(1,(i-1)*35+1:i*35,k)./(Deflator(Country_map(i,2),k)./Deflator(Country_map(i,2),1));
    end
end
sum(CO2)

% CO2_energy=CO2_ENERGY.CO2;

DM=zeros(35*59,2361,9);
for i=1:59
     for j=1:35
         for k=1:9
              DM((i-1)*35+j,:,k)=M((i-1)*35+j,:,k)/(Deflator(Country_map(i,2),k)/Deflator(Country_map(i,2),1));
         end
     end
end

X=DM(1:35*59,2361,:); 
Z=DM(1:35*59,1:35*59,:);
A=zeros(35*59,35*59,9);
h=zeros(35*59,9);
F=DM(1:35*59,35*59+1:35*59+59*5,:);

% F2018=MRIO2018(1:35*59,35*59+1:35*59+59*5);
% Z2018=MRIO2018(1:35*59,1:35*59);

F2=zeros(35*59,59,9);
for i=1:59
   F2(:,i,:)=sum(F(:,(i-1)*5+1:i*5,:),2);
end


%%
%intermediate exports
Export_inter =zeros (59*35,59*35,9);

for i=1:59*35
    for j=1:59*35
        if ceil(i/35)~=ceil(j/35)
            Export_inter (i,j,:) =Z(i,j,:);
        end
    end
end

Export_inter1 =zeros (59*35,59,9);
for i=1:59
      Export_inter1 (:,i,:) =sum(Export_inter(:,(i-1)*35+1:i*35,:),2);
end

%%
for i=1:59 % only service export
 Export_inter1((i-1)*35+1:(i-1)*35+18,:,:)=0;
end

%%final demand 部分      

Export_final= zeros(59*35,59,9);
for i=1:59
   F2(:,i,:)=sum(F(:,(i-1)*5+1:i*5,:),2);
end

for j=1:59*35
   for i=1:59
       if ceil (j/35)~=i
         Export_final(j,i,:)=sum(F(j,(i-1)*5+1:i*5,:),2);
       end
   end
end
for i=1:59 % only service export
 Export_final((i-1)*35+1:(i-1)*35+18,:,:)=0;
 end
Export_T=Export_final+Export_inter1;

save Export_T;
%total trade
Export_trade=zeros(59,59,9);

for i=1:59
    for j=1:59
    Export_trade(i,j,:)=sum(Export_inter1((i-1)*35+1:i*35,j,:))+sum(Export_final((i-1)*35+1:i*35,j,:));
    end
end


%%trade structure
Export_S=zeros(35*59,59,9);
for i=1:59
    for j=1:35
        for k=1:59
              for t=1:9
                 Export_S((i-1)*35+j,k,t)=(Export_final((i-1)*35+j,k,t)+Export_inter1((i-1)*35+j,k,t))/Export_trade(i,k,t);
              end
        end
     end
end

Export_S(isnan(Export_S))=0;
Export_S(isinf(Export_S))=0;
for i=1:59
    for j=1:35
        for k=1:59
              for t=1:9
                checka((i-1)*35+j,k,t)=Export_S((i-1)*35+j,k,t)*Export_trade(i,k,t)-(Export_final((i-1)*35+j,k,t)+Export_inter1((i-1)*35+j,k,t));
              end
        end
     end
end
        
check1=sum(sum(checka(:,:,1)))   % almost 0;


%%


 Arr=zeros(59*35,35,9);

 for i=1:59*35
    for j=1:35
        for k=1:9
         Arr(i,j,k)=Z(i,(ceil(i/35)-1)*35+j,k);
        end
    end
 end

 checkx=sum(sum(Arr(:,:,1)))+sum(sum(F(:,:,1)))+sum(sum(Export_inter1(:,:,1)))-sum(sum(X(:,:,1)))  %  基本为0   -2.9802e-08


 diag_EIO= zeros(59*35,35,9);
%-----------------------------------------------------------output-----------------------------------------------
 Sum_EIO= sum(Z,2)+sum(F,2);

 for i=1:59*35
     for j=1:35
         for k=1:9
            if mod(i-1,35)+1==j
               diag_EIO(i,j,k)=Sum_EIO(i,1,k);
            end
         end
     end
 end
 

 Arr1=zeros(59*35,35,9);

 for i=1:59*35
    for j=1:35
        for k=1:9
         Arr1(i,j,k)= Z(i,(ceil(i/35)-1)*35+j,k)/Sum_EIO((ceil(i/35)-1)*35+j,1,k);
        end
    end
 end
 
Arr1(isnan(Arr1))=0;
Arr1(isinf(Arr1))=0;

B=zeros(59*35,35*59,9);
for K=1:9
        B(:,:,K)= pinv(diag(Sum_EIO(:,:,K)'))*Z(:,:,K);
end
B(isnan(B))=0;
B(isinf(B))=0;
for k=1:9
     LB(:,:,k)=pinv(eye(35*59)-B(:,:,k));
end
LB(isnan(LB))=0;
LB(isinf(LB))=0;


 LL=zeros(59*35,35*59,9);
  for i=1:59
      for k=1:9
         LL(35*(i-1)+1:35*i,35*(i-1)+1:35*i,k)=pinv(diag_EIO(35*(i-1)+1:35*i,1:35,k)-Arr(35*(i-1)+1:35*i,1:35,k));
      end
  end

 L=zeros(59*35,35*59,9);
  for i=1:59
      for k=1:9
         L(35*(i-1)+1:35*i,35*(i-1)+1:35*i,k)=pinv(eye(35)-Arr1(35*(i-1)+1:35*i,1:35,k));
      end
  end

 XL=L(:,:,1)*sum(F2(:,:,1),2);
%%

   for j=1:35*59
        for k=1:9
            h(j,k)=CO2(j,k)./Sum_EIO(j,1,k);%%变为点除
        end
   end

h(isnan(h))=0;
h(isinf(h))=0;

%  
% concordance=xlsread('C:\Users\14718\Desktop\IO\MRIO\Mapping.xlsx','Sheet1');
% %USA 43
% class=concordance(:,2);
% X=V(:,:,1)*LB(:,:,1);
% X2=Sum_EIO(:,:,1);
% X2=X2';
% invertment=zeros(1,8);
% USACO2=zeros(35,8);
% for j=1:8
%     for i=1:59
%         if class(i)==0 && i~=8
%            for k=1:35
%               USACO2(k,j) =  USACO2(k,j)+(V(1,(43-1)*35+1:43*35,j)*LB((43-1)*35+1:43*35,(i-1)*35+k,j)*h((i-1)*35+k,j));
%            end
%         end
%      end
% end
% USA_OD=sum(USACO2(:,:));


% consumption-based emissions
h2=h';
% %only household
% for i=1:59
%    F2(:,i,:)=sum(F(:,(i-1)*5,:));
% end

Flow=zeros(59,9);
 for i=1:9
    Flow(:,i)=h2(i,:)*L(:,:,i)*(F2(:,:,i));
 end
point=[Flow(8,:)
    Flow(11,:)
    Flow(17,:)
    Flow(22,:)
    Flow(25,:)
    Flow(43,:)];

% %%Energy/CO2
CO2_intensity=CO2_energy;
% energy_=zeros(4,35*59,8);
% 
% for k=1:8
%     energy_(:,:,k)=energy(:,:,k)';
% end  
                 
for i=1:35*59
    for k=1:9
       emissionintensity(i,(i-1)*4+1:i*4-1,k)=(CO2_energy(:,i,k)./energy_(1:3,i,k))';
    end
end
emissionintensity(isnan(emissionintensity))=0; 
emissionintensity(isinf(emissionintensity))=0;
CO2_energyi=[sum(CO2_energy(:,(47-1)*35+23,1),2),sum(CO2_energy(:,(47-1)*35+23,5),2),sum(CO2_energy(:,(47-1)*35+23,9),2)];
energy_i=[sum(energy_(1:3,(47-1)*35+23,1),2),sum(energy_(1:3,(47-1)*35+23,5),2),sum(energy_(1:3,(47-1)*35+23,9),2)];

for i=1:35*59
    for k=1:1
       checkco2(:,i,k)=(emissionintensity(i,(i-1)*4+1:i*4-1,k))'.*energy_(1:3,i,k);
    end
end
checkco2=(sum(checkco2))';
checkco22=sum(CO2_energy(:,:,1));

energyshare=energy_./sum(energy_);
for i=1:35*59
    share((i-1)*4+1:i*4,i,:)=energyshare(:,i,:);
end

share(isnan(share))=0;
share(isinf(share))=0;


energytotal=permute(sum(energy_),[2,1,3]);
energyintensity=energytotal./X;
energyintensity(isnan(energyintensity))=0;
energyintensity(isinf(energyintensity))=0;

checka=emissionintensity(:,:,1)*share(:,:,1)*energyintensity(:,:,1);
% c1=emissionintensity(:,:,1);
% c2=share(:,:,1);
% c3=energyintensity(:,:,1).* X;
check0=emissionintensity(:,:,1)*share(:,:,1)*energyintensity(:,:,1).* X;
p0=permute(check0,[1,3,2]);
check2=(checka-h(:,1));  % 验证强度是否是对的 接近0
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
check3=zeros(35*59,59,9);

for i=1:9
  check3(:,:,i)=diag(emissionintensity(:,:,i)*share(:,:,i)*energytotal(:,:,i))*LL(:,:,i)*(Export_inter1(:,:,i)+Export_final(:,:,i));
end
 
sum(sum(check3(:,:,1)))  %  
sum(sum(check3(:,:,2)))-sum(sum(check3(:,:,1)))  %

for i=1:9
  check4(:,:,i)=diag(emissionintensity(:,:,i)*share(:,:,i)*energytotal(:,:,i))*LL(:,:,i)*Export_T(:,:,i);
end

for i=1:59
    for j=1:9
        EXPORT(i,:,j)=sum(check4(35*(i-1)+1:35*i,:,j),1);
    end
end
%save EXPORT;
sum(sum(check4(:,:,2)))-sum(sum(check4(:,:,1)))  % 5.2593e+03

% for i=1:59*35
%     for k=1:9
%       energyintensity(i,1,k)=energytotal(i,1,k)/Sum_EIO(i,1,k);
%     end
% end
% energyintensity(isnan(energyintensity))=0;
% energyintensity(isinf(energyintensity))=0;

% % for i=1:9
% %   check5(:,:,i)=energyintensity(:,:,i).*Sum_EIO(:,:,i)-energytotal(:,:,i);
% % end
% %  sum(sum(check5(:,:,1))) 


for i=1:9
  check6(:,:,i)=diag(emissionintensity(:,:,i)*share(:,:,i)*energyintensity(:,:,i).* Sum_EIO(:,:,i))*LL(:,:,i)*Export_T(:,:,i);
end

  check61=diag(emissionintensity(:,:,1)*share(:,:,1)*energyintensity(:,:,1).* Sum_EIO(:,:,1))*LL(:,:,1);


sum(sum(check6(:,:,2)))-sum(sum(check6(:,:,1)))  % 5.2593e+03

for i=1:9
  check8(:,:,i)=diag(emissionintensity(:,:,i)*share(:,:,i)*energyintensity(:,:,i))*L(:,:,i)*Export_T(:,:,i);
end
 
sum(sum(check8(:,:,2)))-sum(sum(check8(:,:,1))) % 5.2593e+03

Intensity=zeros(35*59,9);
Intensity1=zeros(35*59,9);
CO2_2=CO2';

 for i=1:59
     for k=1:9
       Intensity(35*(i-1)+1:35*i,k)=(CO2_2(k,35*(i-1)+1:35*i)*LL(35*(i-1)+1:35*i,35*(i-1)+1:35*i,k))';    
     end
 end

CO2_2014=zeros(59,35);
for i=1:59
    CO2_2014(i,:)=Intensity((i-1)*35+1:i*35,5)';
end
%    for i=42*35+1:35*43
%      for j=1:9
%         USACO2(i,j) = (Export_final(i,9,j)+ Export_inter1(i,9,j))* Intensity(i,j);
%         %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity1(i,j);
%         %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity(i,j);
%      end
%   end
% USA_CHINA=sum(USACO2(:,:));
% 
%  
%   for i=7*35+1:35*9
%      for j=1:9
%         chianCO2(i,j) = (Export_final(i,43,j)+ Export_inter1(i,43,j))* Intensity(i,j);
%         %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity1(i,j);
%         %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity(i,j);
%      end
%   end
% china_USA=sum(chianCO2(:,:));

 
for i=1:9
  Intensity1(:,i)=sum(diag(emissionintensity(:,:,i)*share(:,:,i)*energyintensity(:,:,i).* Sum_EIO(:,:,i))*LL(:,:,i),2);
end

testCO22017=zeros(35*59,59,9);
 
 for i=1:35*59
     for k=1:59
     for j=1:9
        testCO22017(i,k,j) = (sum(Export_final(i,k,j),2)+ sum(Export_inter1(i,k,j),2))* Intensity(i,j);
        %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity1(i,j);
        %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity(i,j);
     end
     end
 end
 save testCO22017;

testCO22000=zeros(35*59,9);
 
 for i=1:35*59
     for j=1:9
        %testCO22000(i,j) = (sum(Export_final(i,:,j),2)+ sum(Export_inter1(i,:,j),2))* Intensity(i,j);
        testCO22000(i,j) = sum(F2(i,:,j),2)* Intensity1(i,j);
        %testCO22000(i,j) = sum(Export_T(i,:,j),2)* Intensity(i,j);
     end
 end
 
 
 
 
%save testCO22000;
sum(testCO22000(:,:))  %
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%出口
index=[1,3,5,7,9];
emissionintensity=emissionintensity(:,:,index);
share=share(:,:,index);
energyintensity=energyintensity(:,:,index);
L=L(:,:,index);
Export_S=Export_S(:,:,index);
Export_trade=Export_trade(:,:,index);
LL=LL(:,:,index);
Export_T=Export_T(:,:,index);
energytotal=energytotal(:,:,index);
Sum_EIO=Sum_EIO(:,:,index);
CO2=CO2(:,index);

demissionintensity=zeros(35*59,35*59*4,4);
dshare=zeros(4*35*59,35*59,4);
denergyintensity=zeros(35*59,1,4); 
dL=zeros(35*59,35*59,4);
dExport_S=zeros(35*59,59,4);
dExport_trade=zeros(59,59,4);

  for i=1:4
      demissionintensity(:,:,i)=emissionintensity(:,:,i+1)-emissionintensity(:,:,i);
      dshare(:,:,i)=share(:,:,i+1)-share(:,:,i);
      denergyintensity(:,:,i)=energyintensity(:,:,i+1)-energyintensity(:,:,i);
      dL(:,:,i)=L(:,:,i+1)-L(:,:,i);
      dExport_S(:,:,i)=Export_S(:,:,i+1)-Export_S(:,:,i);
      dExport_trade(:,:,i)=Export_trade(:,:,i+1)-Export_trade(:,:,i);
  end
%%
ddemissionintensity1=zeros(35*59,59,4);
ddemissionintensity2=zeros(35*59,59,4);
for t=1:4
    ddemissionintensity1(:,:,t)=diag(demissionintensity(:,:,t)*share(:,:,t)*energytotal(:,:,t))* LL(:,:,t)*Export_T(:,:,t);
    ddemissionintensity2(:,:,t)=diag(demissionintensity(:,:,t)*share(:,:,t+1)*energytotal(:,:,t+1))* LL(:,:,t+1)*Export_T(:,:,t+1);
end
a1=sum(sum(sum(ddemissionintensity1)))  %-104.3851
a2=sum(sum(sum(ddemissionintensity2))) % -113.9599
%%
ddshare1=zeros(35*59,59,4);
ddshare2=zeros(35*59,59,4);
for t=1:4
    ddshare1(:,:,t)=diag(emissionintensity(:,:,t+1)*dshare(:,:,t)*energytotal(:,:,t))* LL(:,:,t)*Export_T(:,:,t);
    ddshare2(:,:,t)=diag(emissionintensity(:,:,t)*dshare(:,:,t)*energytotal(:,:,t+1))* LL(:,:,t+1)*Export_T(:,:,t+1);
end
a3=sum(sum(sum(ddshare1)))  %-104.3851
a4=sum(sum(sum(ddshare2))) % -113.9599
%%
ddenergyintensity1=zeros(35*59,59,4);
ddenergyintensity2=zeros(35*59,59,4);
for t=1:4
    ddenergyintensity1(:,:,t)=diag(emissionintensity(:,:,t+1)*share(:,:,t+1)*denergyintensity(:,:,t).*Sum_EIO(:,:,t))*LL(:,:,t)*Export_T(:,:,t);
    ddenergyintensity2(:,:,t)=diag(emissionintensity(:,:,t)*share(:,:,t)*denergyintensity(:,:,t).*Sum_EIO(:,:,t+1))* LL(:,:,t+1)*Export_T(:,:,t+1);
end
a5=sum(sum(ddenergyintensity1));  %-104.3851
a6=sum(sum(ddenergyintensity2));  % -113.9599
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% production structure
ddL1=zeros(35*59,59,4);
ddL2=zeros(35*59,59,4);

for t=1:4
    ddL1(:,:,t)=diag(emissionintensity(:,:,t+1)*share(:,:,t+1)*energyintensity(:,:,t+1))* dL(:,:,t)*Export_T(:,:,t);
    ddL2(:,:,t)=diag(emissionintensity(:,:,t)*share(:,:,t)*energyintensity(:,:,t))* dL(:,:,t)*Export_T(:,:,t+1);
end 

a8=sum(sum(ddL1));  %-104.3851
a8=sum(sum(ddL2));  % -113.9599
%%
ddExport_S1=zeros(35*59,59,4);
ddExport_S2=zeros(35*59,59,4);

dtemp1=zeros(35*59,59,4);
dtemp2=zeros(35*59,59,4);
for i=1:59
    for j=1:35
        for k=1:59
              for t=1:4
                 dtemp1((i-1)*35+j,k,t)=dExport_S((i-1)*35+j,k,t)*Export_trade(i,k,t);
                 dtemp2((i-1)*35+j,k,t)=dExport_S((i-1)*35+j,k,t)*Export_trade(i,k,t+1);
              end
        end
     end
end
for t=1:4
    ddExport_S1(:,:,t)=diag(CO2(:,t+1))* LL(:,:,t+1)*dtemp1(:,:,t);
    ddExport_S2(:,:,t)=diag(CO2(:,t))* LL(:,:,t)*dtemp2(:,:,t);
end 
a9=sum(sum(ddExport_S1));  %-104.3851
a10=sum(sum(ddExport_S2)); % -113.9599
%%
ddExport_trade1=zeros(35*59,59,4);
ddExport_trade2=zeros(35*59,59,4);

dtemp1=zeros(35*59,59,4);
dtemp2=zeros(35*59,59,4);
for i=1:59
    for j=1:35
        for k=1:59
              for t=1:4
                 dtemp1((i-1)*35+j,k,t)=Export_S((i-1)*35+j,k,t+1)*dExport_trade(i,k,t);
                 dtemp2((i-1)*35+j,k,t)=Export_S((i-1)*35+j,k,t)*dExport_trade(i,k,t);
              end
        end
     end
end

for t=1:4
    ddExport_trade1(:,:,t)=diag(CO2(:,t+1))* LL(:,:,t+1)*dtemp1(:,:,t);
    ddExport_trade2(:,:,t)=diag(CO2(:,t))* LL(:,:,t)*dtemp2(:,:,t);
end 
sum(sum(ddExport_trade1));  %-104.3851
sum(sum(ddExport_trade2));  % -113.9599

 ddemissionintensity=(ddemissionintensity1+ddemissionintensity2)/2;
 ddshare=(ddshare1+ddshare2)/2;
 ddenergyintensity=(ddenergyintensity1+ddenergyintensity2)/2;
 ddL=(ddL1+ddL2)/2;
 ddExport_S=(ddExport_S1+ddExport_S2)/2;
 ddExport_trade=(ddExport_trade1+ddExport_trade2)/2;
 check8=sum(sum(ddemissionintensity))+sum(sum(ddshare))+sum(sum(ddenergyintensity))+sum(sum(ddL))+sum(sum(ddExport_S))+sum(sum(ddExport_trade));

 id=[1,3,5,7,9];
 check4=check4(:,:,id);
for i=1:4
  check9(:,:,i)=check4(:,:,i+1)-check4(:,:,i);
end
 SDA_emissionintensity=zeros(59,59,4);
 SDA_share=zeros(59,59,4);
 SDA_energyintensity=zeros(59,59,4);
 SDA_L=zeros(59,59,4);
 SDA_Export_S=zeros(59,59,4);
 SDA_Export_trade=zeros(59,59,4);
 total=zeros(59,59,4);
 for i=1:59
    SDA_emissionintensity(i,:,:)=sum(ddemissionintensity((i-1)*35+1:i*35,:,:)); 
    SDA_share(i,:,:)=sum(ddshare((i-1)*35+1:i*35,:,:)); 
    SDA_energyintensity(i,:,:)=sum(ddenergyintensity((i-1)*35+1:i*35,:,:)); 
    SDA_L(i,:,:)=sum(ddL((i-1)*35+1:i*35,:,:)); 
    SDA_Export_S(i,:,:)=sum(ddExport_S((i-1)*35+1:i*35,:,:)); 
    SDA_Export_trade(i,:,:)=sum(ddExport_trade((i-1)*35+1:i*35,:,:)); 
    total(i,:,:)=sum(check9((i-1)*35+1:i*35,:,:)); 
 end

 %分解的判断check10
 check10=zeros(59,59,4);
 for k=1:4
    T=SDA_emissionintensity(:,:,k)+SDA_share(:,:,k)+SDA_energyintensity(:,:,k)+SDA_L(:,:,k)+SDA_Export_S(:,:,k)+ SDA_Export_trade(:,:,k); 
    check10(:,:,k)=total(:,:,k)-T;
 end
 
 % SDA2016_2017_D=[Japan_SDA_D(:,3),Singapore_SDA_D(:,3),Germany_SDA_D(:,3),INDIA_SDA_D(:,3),France_SDA_D(:,3),UK_SDA_D(:,3),Australia_SDA_D(:,3),Canada_SDA_D(:,3),USA_SDA_D(:,3),Italy_SDA_D(:,3),China_SDA_D(:,3)];

 %% S-N
 Mapping_country=xlsread('C:\Users\14718\Desktop\IO\MRIO\Mapping.xlsx','MAPPING_EU28');
 EU_index=Mapping_country(:,1);% USA 43 CHINA 4 INDIA 22
 developing_index=[8;22;54;56];
 developed_index=[1;25;43];
 %%
%DD DP PD PP
concordance=xlsread('C:\Users\14718\Desktop\IO\MRIO\Mapping.xlsx','Sheet1');
class=concordance(:,2);
%developed to developed
DD_SDA=[permute(sum(sum(SDA_emissionintensity(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_share(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_energyintensity(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_L(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_Export_S(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_Export_trade(class==1,class==1,:))),[1,3,2])
        permute(sum(sum(total(class==1,class==1,:))),[1,3,2])];

%developed to developing
DP_SDA=[permute(sum(sum(SDA_emissionintensity(class==1,class==0,:))),[1,3,2])
        permute(sum(sum(SDA_share(class==1,class==0,:))),[1,3,2])
        permute(sum(sum(SDA_energyintensity(class==1,class==0,:))),[1,3,2])
        permute(sum(sum(SDA_L(class==1,class==0,:))),[1,3,2])
        permute(sum(sum(SDA_Export_S(class==1,class==0,:))),[1,3,2])
        permute(sum(sum(SDA_Export_trade(class==1,class==0,:))),[1,3,2]) 
        permute(sum(sum(total(class==1,class==0,:))),[1,3,2])];
%developing to developed
PD_SDA=[permute(sum(sum(SDA_emissionintensity(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_share(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_energyintensity(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_L(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_Export_S(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(SDA_Export_trade(class==0,class==1,:))),[1,3,2])
        permute(sum(sum(total(class==0,class==1,:))),[1,3,2])];
%developing to developing
PP_SDA=[permute(sum(sum(SDA_emissionintensity(class==0,class==0,:))),[1,3,2]) % emission intensity
        permute(sum(sum(SDA_share(class==0,class==0,:))),[1,3,2]) % energy structure
        permute(sum(sum(SDA_energyintensity(class==0,class==0,:))),[1,3,2]) % energy intensity
        permute(sum(sum(SDA_L(class==0,class==0,:))),[1,3,2])  % production structure
        permute(sum(sum(SDA_Export_S(class==0,class==0,:))),[1,3,2]) % trade structure
        permute(sum(sum(SDA_Export_trade(class==0,class==0,:))),[1,3,2]) % trade volume
        permute(sum(sum(total(class==0,class==0,:))),[1,3,2])]; % total emission chage    
    
VC_SDA=[permute(sum(SDA_emissionintensity(48,8,:),2),[1,3,2])+permute(sum(SDA_share(48,8,:),2),[1,3,2])+permute(sum(SDA_energyintensity(48,8,:),2),[1,3,2])
        permute(sum(SDA_L(48,8,:),2),[1,3,2])
        permute(sum(SDA_Export_S(48,8,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(48,8,:),2),[1,3,2])
        permute(sum(total(48,8,:),2),[1,3,2])];
    
IC_SDA=[permute(sum(SDA_emissionintensity(22,8,:),2),[1,3,2])+permute(sum(SDA_share(22,8,:),2),[1,3,2])+permute(sum(SDA_energyintensity(22,8,:),2),[1,3,2])
        permute(sum(SDA_L(22,8,:),2),[1,3,2])
        permute(sum(SDA_Export_S(22,8,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(22,8,:),2),[1,3,2])
        permute(sum(total(22,8,:),2),[1,3,2])];
    
TI_SDA=[permute(sum(SDA_emissionintensity(47,21,:),2),[1,3,2])+permute(sum(SDA_share(47,21,:),2),[1,3,2])+permute(sum(SDA_energyintensity(47,21,:),2),[1,3,2])
        permute(sum(SDA_L(47,21,:),2),[1,3,2])
        permute(sum(SDA_Export_S(47,21,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(47,21,:),2),[1,3,2])
        permute(sum(total(47,21,:),2),[1,3,2])];
    
VJ_SDA=[permute(sum(SDA_emissionintensity(48,25,:),2),[1,3,2])+permute(sum(SDA_share(48,25,:),2),[1,3,2])+permute(sum(SDA_energyintensity(48,25,:),2),[1,3,2])
        permute(sum(SDA_L(48,25,:),2),[1,3,2])
        permute(sum(SDA_Export_S(48,25,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(48,25,:),2),[1,3,2])
        permute(sum(total(48,25,:),2),[1,3,2])];
    
CU_SDA=[permute(sum(SDA_emissionintensity(8,43,:),2),[1,3,2])+permute(sum(SDA_share(8,43,:),2),[1,3,2])+permute(sum(SDA_energyintensity(8,43,:),2),[1,3,2])
        permute(sum(SDA_L(8,43,:),2),[1,3,2])
        permute(sum(SDA_Export_S(8,43,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(8,43,:),2),[1,3,2])
        permute(sum(total(8,43,:),2),[1,3,2])];
    
UC_SDA=[permute(sum(SDA_emissionintensity(43,8,:),2),[1,3,2])+permute(sum(SDA_share(43,8,:),2),[1,3,2])+permute(sum(SDA_energyintensity(43,8,:),2),[1,3,2])
        permute(sum(SDA_L(43,8,:),2),[1,3,2])
        permute(sum(SDA_Export_S(43,8,:),2),[1,3,2])
        permute(sum(SDA_Export_trade(43,8,:),2),[1,3,2])
        permute(sum(total(43,8,:),2),[1,3,2])];