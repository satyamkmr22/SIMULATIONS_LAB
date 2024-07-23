clc;
clear all;
data1=readmatrix("data1.xlsx");
data2=readmatrix('data_caliberation.xlsx');
caliberation_absorbance=data2(:,1);
caliberation_distilled_water=data2(:,2);
caliberation_KMnO4=data2(:,3);
weight_KMnO4=5;
distilled_water=10;
M_mass_KMnO4=158.034;
Concentration_KMnO4=weight_KMnO4/(M_mass_KMnO4*distilled_water*1000);
final_conc_caliberation=Concentration_KMnO4.*caliberation_KMnO4./(caliberation_distilled_water+caliberation_KMnO4);
time=data1(:,1);
absorbance_low=data1(:,3);
absorbance_mid=data1(:,4);
absorbance_top=data1(:,2);
figure(1)
plot(time,absorbance_low);
hold on;
plot(time,absorbance_mid);
hold on;
plot(time,absorbance_top);
xlabel('time');
ylabel('absorbance');
legend('lowest point','middle point','top');
figure(2)
plot(final_conc_caliberation,caliberation_absorbance);
xlabel('concentration KMnO4');
ylabel('absorbance');
