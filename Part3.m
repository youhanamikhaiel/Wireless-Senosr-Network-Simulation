clear all
clc
n = 100; %number of sensor nodes
l = 100; %length of area
w = 100; %width of area
Ei = ones(1,100)*2*(10^9);  %Initial energy for single node in nJ
eps_short = 10;
eps_long  = 0.0013;
Eelec = 50; % in nJ/bit
Eagg = 50; % in nJ/bit
Data = 500*8; % in bits
Overhead = 125*8; % in bits
kbits = Data + Overhead; % Total number of bits sent in one frame
kbitsagg = 500*8;
d0 = sqrt(eps_short/eps_long); % Critical distance
xy = randi([0 n],2,n); % Nodes' locations
sink = [w/2 l/2]; % Sink location
NoCHs = 0.05*n;
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l])
xlim([0 w])
hold on
scatter(sink(1),sink(2),'filled')
title('Original network topology')
dist = [];
for i=1:n
    dist(i) = sqrt(((xy(1,i)-sink(1))^2) + ((xy(2,i)-sink(2))^2));
end
hhg = 1;
for i=1:n
    Ei1(i) = Ei(i) + (hhg*5);
    hhg = hhg + 1;
end
[no_cl, cen] = kmeans(xy',NoCHs);
cen = cen';
cluster1 = xy(:,find(no_cl==1));
cluster2 = xy(:,find(no_cl==2));
cluster3 = xy(:,find(no_cl==3));
cluster4 = xy(:,find(no_cl==4));
cluster5 = xy(:,find(no_cl==5));
figure
hold on
scatter(cluster1(1,:),cluster1(2,:),'filled')
scatter(cluster2(1,:),cluster2(2,:),'filled')
scatter(cluster3(1,:),cluster3(2,:),'filled')
scatter(cluster4(1,:),cluster4(2,:),'filled')
scatter(cluster5(1,:),cluster5(2,:),'filled')
title('Network topology after clustering')
C = 5; % Number of cycles before electing new CH
notdead1 = ones(1,n);
for j=1:1000
%     for k=1:length(cluster1)
%         cen_dist1(k) = sqrt(((cluster1(1,k)-cen(1,1))^2) + ((cluster1(2,k)-cen(2,1))^2));
%     end
%     CH1 = find(cen_dist==min(cen_dist));
    CH1 = find(Ei1==max(Ei1(no_cl==1)));
    CH2 = find(Ei1==max(Ei1(no_cl==2)));
    CH3 = find(Ei1==max(Ei1(no_cl==3)));
    CH4 = find(Ei1==max(Ei1(no_cl==4)));
    CH5 = find(Ei1==max(Ei1(no_cl==5)));
    CHeads = [CH1, CH2, CH3, CH4, CH5];
    len_clusters = [length(cluster1), length(cluster2), length(cluster3), length(cluster4), length(cluster5)];
    for i=1:n
            dist_chs_nodes(i) = sqrt(((xy(1,i)-xy(1,CHeads(no_cl(i))))^2) + ((xy(2,i)-xy(1,CHeads(no_cl(i))))^2));
            CHs_dist_sink(i) = dist(CHeads(no_cl(i)));
    end
    for i=1:C
        for e=1:n
            energy_loss1 = kbits*Eelec + kbits*eps_short*(dist_chs_nodes(e)^2);
            energy_loss2 = kbits*Eelec + kbits*eps_short*(CHs_dist_sink(e)^2) + kbits*Eagg;
            if Ei1(CHeads(no_cl(e))) > C*energy_loss2
                if Ei1(e) > energy_loss1
                    Ei1(e) = Ei1(e) - energy_loss1;
                    Ei1(CHeads(no_cl(e))) = Ei1(CHeads(no_cl(e))) - energy_loss2;
                else
                    if notdead1(e) == 1
                        notdead1(e) = 0;
                        death1(e) = j*C;
                    end
                    
                end
            else
                if notdead1(CHeads(no_cl(e))) == 1
                    notdead1(CHeads(no_cl(e))) = 0;
                    death1(CHeads(no_cl(e))) = j*C;
                end
            end
        end
    end
end
death1(death1==0)=15;
death1(death1==5)=15;
death1(death1==10)=20;
T111 = min(death1);
notedead3 = [];
for j=1:max(death1)
    nn = 00;
    for i=1:n
        if death1(i) >= j
            nn = nn + 1;
        else
            nn = nn;
        end
    end
    notdead3(j) = nn;
end
notdead3 = [notdead3 0];
no_cycles = 1:max(death1)+1;

figure
plot(no_cycles, notdead3)
xlabel('Number of cycles')
ylabel('Number of living nodes')
ylim([0 105])
hold on 
plot(ones(n+1)*T111,0:n)
title('Number of active nodes vs number of cycles at C=5')

%%

Ei1=Ei;
notdead1 = ones(1,n);
for j=1:T111/5
    %     for k=1:length(cluster1)
    %         cen_dist1(k) = sqrt(((cluster1(1,k)-cen(1,1))^2) + ((cluster1(2,k)-cen(2,1))^2));
    %     end
    %     CH1 = find(cen_dist==min(cen_dist));
    CH1 = find(Ei1==max(Ei1(no_cl==1)));
    CH2 = find(Ei1==max(Ei1(no_cl==2)));
    CH3 = find(Ei1==max(Ei1(no_cl==3)));
    CH4 = find(Ei1==max(Ei1(no_cl==4)));
    CH5 = find(Ei1==max(Ei1(no_cl==5)));
    CHeads = [CH1, CH2, CH3, CH4, CH5];
    len_clusters = [length(cluster1), length(cluster2), length(cluster3), length(cluster4), length(cluster5)];
    for i=1:n
        dist_chs_nodes(i) = sqrt(((xy(1,i)-xy(1,CHeads(no_cl(i))))^2) + ((xy(2,i)-xy(1,CHeads(no_cl(i))))^2));
        CHs_dist_sink(i) = dist(CHeads(no_cl(i)));
    end
    for i=1:C
        for e=1:n
            energy_loss1 = kbits*Eelec + kbits*eps_short*(dist_chs_nodes(e)^2);
            energy_loss2 = kbits*Eelec + kbits*eps_short*(CHs_dist_sink(e)^2) + kbits*Eagg;
            if Ei1(CHeads(no_cl(e))) > C*energy_loss2
                if Ei1(e) > energy_loss1
                    Ei1(e) = Ei1(e) - energy_loss1;
                    Ei1(CHeads(no_cl(e))) = Ei1(CHeads(no_cl(e))) - energy_loss2;
                else
                    if notdead1(e) == 1
                        notdead1(e) = 0;
                        death1(e) = j*C;
                    end
                    
                end
            else
                if notdead1(CHeads(no_cl(e))) == 1
                    notdead1(CHeads(no_cl(e))) = 0;
                    death1(CHeads(no_cl(e))) = j*C;
                end
            end
        end
    end
end
figure
stem(1:n,Ei1)
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T1 cycle in case 1')


%%
yyt = 1;
for C=1:15
    Ei1=Ei;
    notdead1 = ones(1,n);
    for j=1:1000
        %     for k=1:length(cluster1)
        %         cen_dist1(k) = sqrt(((cluster1(1,k)-cen(1,1))^2) + ((cluster1(2,k)-cen(2,1))^2));
        %     end
        %     CH1 = find(cen_dist==min(cen_dist));
        CH1 = find(Ei1==max(Ei1(no_cl==1)));
        CH2 = find(Ei1==max(Ei1(no_cl==2)));
        CH3 = find(Ei1==max(Ei1(no_cl==3)));
        CH4 = find(Ei1==max(Ei1(no_cl==4)));
        CH5 = find(Ei1==max(Ei1(no_cl==5)));
        CHeads = [CH1, CH2, CH3, CH4, CH5];
        len_clusters = [length(cluster1), length(cluster2), length(cluster3), length(cluster4), length(cluster5)];
        for i=1:n
            dist_chs_nodes(i) = sqrt(((xy(1,i)-xy(1,CHeads(no_cl(i))))^2) + ((xy(2,i)-xy(1,CHeads(no_cl(i))))^2));
            CHs_dist_sink(i) = dist(CHeads(no_cl(i)));
        end
        for i=1:C
            for e=1:n
                energy_loss1 = kbits*Eelec + kbits*eps_short*(dist_chs_nodes(e)^2);
                energy_loss2 = kbits*Eelec + kbits*eps_short*(CHs_dist_sink(e)^2) + kbits*Eagg;
                if Ei1(CHeads(no_cl(e))) > C*energy_loss2
                    if Ei1(e) > energy_loss1
                        Ei1(e) = Ei1(e) - energy_loss1;
                        Ei1(CHeads(no_cl(e))) = Ei1(CHeads(no_cl(e))) - energy_loss2;
                    else
                        if notdead1(e) == 1
                            notdead1(e) = 0;
                            death1(e) = j*C;
                        end
                        
                    end
                else
                    if notdead1(CHeads(no_cl(e))) == 1
                        notdead1(CHeads(no_cl(e))) = 0;
                        death1(CHeads(no_cl(e))) = j*C;
                    end
                end
            end
        end
    end
    T111array(yyt) = min(death1);
    yyt = yyt + 1;
end
Carray = [1:15];
figure
plot(Carray,T111array)
title('Changing number of cycles versus death of first node T1')
ylim([5 20])
xlabel('Number of cycles before new cluster head election')
ylabel('Death of first node T1')



%% Part E
Eii1 = Ei;
internodes = [];
no_internodes = 5;
Einter = (4*(10^9))*ones(1,no_internodes); % Initial energy of intermediate nodes in nJ
R = 25;
for i=1:no_internodes
    internodes(1,i)= 50 + R*sind(i*(360/no_internodes));
    internodes(2,i)= 50 + R*cosd(i*(360/no_internodes));
end
figure
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l])
xlim([0 w])
hold on
scatter(sink(1),sink(2),'filled')
hold on
scatter(internodes(1,:),internodes(2,:), 'filled')
angle = 0:0.01:2*pi;
xc = R*cos(angle)+sink(1);
yc = R*sin(angle)+sink(2);
hold on
plot(xc,yc)
title('Network realization for 5 uniformly distributed CHs on a circle, R = 25 m')
internodes_index = [];
dist_internodes = [];
for i=1:n
    for j=1:no_internodes
        dist_internodes(j,i) = sqrt(((xy(1,i)-internodes(1,j))^2) + ((xy(2,i)-internodes(2,j))^2));
    end
end
min_dist_internodes = min(dist_internodes);
for i=1:n
    internodes_index(i) = find(dist_internodes(:,i)== min_dist_internodes(i));
end
death_inter = [];
not_dead2 = ones(1,n);
energy_loss_inter = kbits*Eelec + kbits*eps_short*(R^2) + 2*kbits*Eagg;
for cycles=1:5000
    for i=1:n
        energy_loss_short_inter = kbits*Eelec + kbits*eps_short*(min_dist_internodes(i)^2);
        energy_loss_long_inter = kbits*Eelec + kbits*eps_long*(min_dist_internodes(i)^4);
        energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2);
        energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4);
        if Einter(internodes_index(i)) > energy_loss_inter
            if min_dist_internodes < d0
                if Eii1(i) > energy_loss_short_inter
                    Eii1(i) = Eii1(i) - energy_loss_short_inter;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            else
                if Eii1(i) > energy_loss_long_inter
                    Eii1(i) = Eii1(i) - energy_loss_long_inter;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            end
            Einter(internodes_index(i)) = Einter(internodes_index(i)) - energy_loss_inter;
        else
            if min_dist_internodes < d0
                if Eii1(i) > energy_loss_short
                    Eii1(i) = Eii1(i) - energy_loss_short;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            else
                if Eii1(i) > energy_loss_long
                    Eii1(i) = Eii1(i) - energy_loss_long;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            end
        end
    end
end

T11 = min(death_inter);
notedead2 = [];
for j=1:max(death_inter)
    nn = 00;
    for i=1:n
        if death_inter(i) >= j
            nn = nn + 1;
        else
            nn = nn;
        end
    end
    notdead2(j) = nn;
end
notdead2 = [notdead2 0];
no_cycles = 1:max(death_inter)+1;

figure
plot(no_cycles, notdead2)
xlabel('Number of cycles')
ylabel('Number of living nodes')
ylim([0 105])
hold on 
plot(ones(n+1)*T11,0:n)
title('Number of active nodes vs number of cycles at R=25')

%%
Eii1 = Ei;
internodes = [];
no_internodes = 5;
Einter = (4*(10^9))*ones(1,no_internodes); % Initial energy of intermediate nodes in nJ
R = 25;
for i=1:no_internodes
    internodes(1,i)= 50 + R*sind(i*(360/no_internodes));
    internodes(2,i)= 50 + R*cosd(i*(360/no_internodes));
end
% figure
% scatter(xy(1,:),xy(2,:),'filled')
% ylim([0 l])
% xlim([0 w])
% hold on
% scatter(sink(1),sink(2),'filled')
% hold on
% scatter(internodes(1,:),internodes(2,:), 'filled')
% angle = 0:0.01:2*pi;
% xc = R*cos(angle)+sink(1);
% yc = R*sin(angle)+sink(2);
% hold on
% plot(xc,yc)
% title('Network realization for 5 uniformly distributed CHs on a circle, R = 25 m')
internodes_index = [];
dist_internodes = [];
for i=1:n
    for j=1:no_internodes
        dist_internodes(j,i) = sqrt(((xy(1,i)-internodes(1,j))^2) + ((xy(2,i)-internodes(2,j))^2));
    end
end
min_dist_internodes = min(dist_internodes);
for i=1:n
    internodes_index(i) = find(dist_internodes(:,i)== min_dist_internodes(i));
end
death_inter = [];
not_dead2 = ones(1,n);
energy_loss_inter = kbits*Eelec + kbits*eps_short*(R^2) + 2*kbits*Eagg;
for cycles=1:T11
    for i=1:n
        energy_loss_short_inter = kbits*Eelec + kbits*eps_short*(min_dist_internodes(i)^2);
        energy_loss_long_inter = kbits*Eelec + kbits*eps_long*(min_dist_internodes(i)^4);
        energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2);
        energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4);
        if Einter(internodes_index(i)) > energy_loss_inter
            if min_dist_internodes < d0
                if Eii1(i) > energy_loss_short_inter
                    Eii1(i) = Eii1(i) - energy_loss_short_inter;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            else
                if Eii1(i) > energy_loss_long_inter
                    Eii1(i) = Eii1(i) - energy_loss_long_inter;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            end
            Einter(internodes_index(i)) = Einter(internodes_index(i)) - energy_loss_inter;
        else
            if min_dist_internodes < d0
                if Eii1(i) > energy_loss_short
                    Eii1(i) = Eii1(i) - energy_loss_short;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            else
                if Eii1(i) > energy_loss_long
                    Eii1(i) = Eii1(i) - energy_loss_long;
                else
                    if not_dead2(i) == 1
                        not_dead2(i) = 0;
                        death_inter(i) = cycles;
                    end
                end
            end
        end
    end
end
figure
stem(Eii1)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T1 cycles for case 2')

%% Part F
nhn = 1;
for R=5:5:80
    Eii1 = Ei;
    internodes = [];
    Einter = (4*(10^9))*ones(1,no_internodes); % Initial energy of intermediate nodes in nJ
    no_internodes = 5;
    for i=1:no_internodes
        internodes(1,i)= 50 + R*sind(i*(360/no_internodes));
        internodes(2,i)= 50 + R*cosd(i*(360/no_internodes));
    end
    % figure
    % scatter(xy(1,:),xy(2,:),'filled')
    % ylim([0 l])
    % xlim([0 w])
    % hold on
    % scatter(sink(1),sink(2),'filled')
    % hold on
    % scatter(internodes(1,:),internodes(2,:), 'filled')
    % angle = 0:0.01:2*pi;
    % xc = R*cos(angle)+sink(1);
    % yc = R*sin(angle)+sink(2);
    % hold on
    % plot(xc,yc)
    % title('Network realization for 5 uniformly distributed CHs on a circle, R = 25 m')
    internodes_index = [];
    dist_internodes = [];
    for i=1:n
        for j=1:no_internodes
            dist_internodes(j,i) = sqrt(((xy(1,i)-internodes(1,j))^2) + ((xy(2,i)-internodes(2,j))^2));
        end
    end
    min_dist_internodes = min(dist_internodes);
    for i=1:n
        internodes_index(i) = find(dist_internodes(:,i)== min_dist_internodes(i));
    end
    death_inter = [];
    not_dead2 = ones(1,n);
    energy_loss_inter = kbits*Eelec + kbits*eps_short*(R^2) + 2*kbits*Eagg;
    for cycles=1:5000
        for i=1:n
            energy_loss_short_inter = kbits*Eelec + kbits*eps_short*(min_dist_internodes(i)^2);
            energy_loss_long_inter = kbits*Eelec + kbits*eps_long*(min_dist_internodes(i)^4);
            energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2);
            energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4);
            if Einter(internodes_index(i)) > energy_loss_inter
                if min_dist_internodes < d0
                    if Eii1(i) > energy_loss_short_inter
                        Eii1(i) = Eii1(i) - energy_loss_short_inter;
                    else
                        if not_dead2(i) == 1
                            not_dead2(i) = 0;
                            death_inter(i) = cycles;
                        end
                    end
                else
                    if Eii1(i) > energy_loss_long_inter
                        Eii1(i) = Eii1(i) - energy_loss_long_inter;
                    else
                        if not_dead2(i) == 1
                            not_dead2(i) = 0;
                            death_inter(i) = cycles;
                        end
                    end
                end
                Einter(internodes_index(i)) = Einter(internodes_index(i)) - energy_loss_inter;
            else
                if min_dist_internodes < d0
                    if Eii1(i) > energy_loss_short
                        Eii1(i) = Eii1(i) - energy_loss_short;
                    else
                        if not_dead2(i) == 1
                            not_dead2(i) = 0;
                            death_inter(i) = cycles;
                        end
                    end
                else
                    if Eii1(i) > energy_loss_long
                        Eii1(i) = Eii1(i) - energy_loss_long;
                    else
                        if not_dead2(i) == 1
                            not_dead2(i) = 0;
                            death_inter(i) = cycles;
                        end
                    end
                end
            end
        end
    end
    
    T11array(nhn) = min(death_inter);
    T33array(nhn) = max(death_inter);
    nhn = nhn + 1;
end
figure
Rloop = [5:5:80];
subplot(2,1,1)
plot(Rloop, T11array)
xlabel('Radius of the circle')
ylabel('Number of cycles')
title('T1 versus changing the radius')
ylim([0 20])
subplot(2,1,2)
plot(Rloop, T33array)
xlabel('Radius of the circle')
ylabel('Number of cycles')
title('T3 versus changing the radius')







