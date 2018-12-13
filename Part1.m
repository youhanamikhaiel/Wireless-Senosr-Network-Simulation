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
k = Data + Overhead;
d0 = sqrt(eps_short/eps_long);
xy = randi([0 n],2,n);
sink = [w/2 l/2];
sink2 = [0 0];
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l])
xlim([0 w])
hold on
scatter(sink(1),sink(2),'filled')
dist = [];
for i=1:n
    dist(i) = sqrt(((xy(1,i)-sink(1))^2) + ((xy(2,i)-sink(2))^2));
end
Ei1 = Ei;
death = [];
for i=1:n
    m = 0;
    energy_loss_short = k*Eelec + k*eps_short*(dist(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist(i)^4) + k*Eagg;
    if dist(i) < d0
        while energy_loss_short < Ei1(i)
            Ei1(i) = Ei1(i) - energy_loss_short;
            m = m + 1;
        end
        death(i) = m;
    else
        while energy_loss_long < Ei1(i)
            Ei1(i) = Ei1(i) - energy_loss_long;
            m = m + 1;
        end
        death(i) = m;
    end
end
%number of cycles before the death of the first node
T1 = min(death);
rt11 = find(death==min(death));
%number of cycles before the death of the last node
T3 = max(death);
for i=1:length(rt11)
    hold on
    scatter(xy(1,rt11(i)),xy(2,rt11(i)),'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor',[0.5 0 0],'LineWidth',3)
    title('Generated network topology')
end
notdead = [];
for j=1:max(death)
    nn = 00;
    for i=1:n
        if death(i) >= j
            nn = nn + 1;
        else
            nn = nn;
        end
    end
    notdead(j) = nn;
end
notdead2 = [notdead 0];
cycles = 1:max(death)+1;
figure
plot(cycles, notdead2)
xlabel('Number of cycles')
ylabel('Number of living nodes')
ylim([0 105])

kk=1;
while notdead(kk)>50
    kk = kk + 1;
end
%number of cycles before the death of half the nodes
T2 = kk;
hold on 
plot(ones(n+1)*T1,0:n)
hold on 
plot(ones(n+1)*T2,0:n)
hold on 
plot(ones(n+1)*T3,0:n)

%Plot the remaining energies after T1 cycle
Ei2 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist(i)^4) + k*Eagg;
    if dist(i) < d0
        for ee=1:T1
            Ei2(i) = Ei2(i) - energy_loss_short;
        end
    else
        for ee=1:T1
            Ei2(i) = Ei2(i) - energy_loss_long; 
        end
    end
end
figure
stem(1:n,Ei2)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T1 cycle')

%Plot the remaining energies after T2 cycles
Ei3 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist(i)^4) + k*Eagg;
    if dist(i) < d0
        for j=1:T2
            if energy_loss_short <= Ei3(i)
                Ei3(i) = Ei3(i) - energy_loss_short;
            end
        end
    else
        for j=1:T2
            if energy_loss_long <= Ei3(i)
                Ei3(i) = Ei3(i) - energy_loss_long;
            end
        end
    end
end
figure
stem(1:n,Ei3)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T2 cycle')

%Plot the remaining energies after T2 cycles
Ei4 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist(i)^4) + k*Eagg;
    if dist(i) < d0
        for j=1:T3
            if energy_loss_short <= Ei4(i)
                Ei4(i) = Ei4(i) - energy_loss_short;
            end
        end
    else
        for j=1:T3
            if energy_loss_long <= Ei4(i)
                Ei4(i) = Ei4(i) - energy_loss_long;
            end
        end
    end
end
figure
stem(1:n,Ei4)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title(['Remaining energies after T3 cycles, T3=' num2str(T3) ' cycles'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l])
xlim([0 w])
hold on
scatter(sink2(1),sink2(2),'filled')
dist = [];
for i=1:n
    dist2(i) = sqrt(((xy(1,i)-sink2(1))^2) + ((xy(2,i)-sink2(2))^2));
end
Ei1 = Ei;
death = [];
for i=1:n
    m = 0;
    energy_loss_short = k*Eelec + k*eps_short*(dist2(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist2(i)^4) + k*Eagg;
    if dist2(i) < d0
        while energy_loss_short < Ei1(i)
            Ei1(i) = Ei1(i) - energy_loss_short;
            m = m + 1;
        end
        death(i) = m;
    else
        while energy_loss_long < Ei1(i)
            Ei1(i) = Ei1(i) - energy_loss_long;
            m = m + 1;
        end
        death(i) = m;
    end
end
%number of cycles before the death of the first node
T1 = min(death);
rt11 = find(death==min(death));
%number of cycles before the death of the last node
T3 = max(death);
for i=1:length(rt11)
    hold on
    scatter(xy(1,rt11(i)),xy(2,rt11(i)),'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor',[0.5 0 0],'LineWidth',3)
    title('Generated network topology - corner sink node')
end
notdead = [];
for j=1:max(death)
    nn = 00;
    for i=1:n
        if death(i) >= j
            nn = nn + 1;
        else
            nn = nn;
        end
    end
    notdead(j) = nn;
end
notdead2 = [notdead 0];
cycles = 1:max(death)+1;
figure
plot(cycles, notdead2)
xlabel('Number of cycles')
ylabel('Number of living nodes')
ylim([0 105])
title('For corner sink node')
kk=1;
while notdead(kk)>50
    kk = kk + 1;
end
%number of cycles before the death of half the nodes
T2 = kk;
hold on 
plot(ones(n+1)*T1,0:n)
hold on 
plot(ones(n+1)*T2,0:n)
hold on 
plot(ones(n+1)*T3,0:n)


%Plot the remaining energies after T1 cycle
Ei2 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist2(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist2(i)^4) + k*Eagg;
    if dist2(i) < d0
        for ee=1:T1
            Ei2(i) = Ei2(i) - energy_loss_short;
        end
    else
        for ee=1:T1
            Ei2(i) = Ei2(i) - energy_loss_long; 
        end
    end
end
figure
stem(1:n,Ei2)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T1 cycle - corner sink node')

%Plot the remaining energies after T2 cycles
Ei3 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist2(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist2(i)^4) + k*Eagg;
    if dist2(i) < d0
        for j=1:T2
            if energy_loss_short <= Ei3(i)
                Ei3(i) = Ei3(i) - energy_loss_short;
            end
        end
    else
        for j=1:T2
            if energy_loss_long <= Ei3(i)
                Ei3(i) = Ei3(i) - energy_loss_long;
            end
        end
    end
end
figure
stem(1:n,Ei3)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title('Remaining energies after T2 cycle - corner sink node')

%Plot the remaining energies after T2 cycles
Ei4 = Ei;
for i=1:n
    energy_loss_short = k*Eelec + k*eps_short*(dist2(i)^2) + k*Eagg;
    energy_loss_long = k*Eelec + k*eps_long*(dist2(i)^4) + k*Eagg;
    if dist2(i) < d0
        for j=1:T3
            if energy_loss_short <= Ei4(i)
                Ei4(i) = Ei4(i) - energy_loss_short;
            end
        end
    else
        for j=1:T3
            if energy_loss_long <= Ei4(i)
                Ei4(i) = Ei4(i) - energy_loss_long;
            end
        end
    end
end
figure
stem(1:n,Ei4)
hold on 
plot(1:n,Ei)
ylim([0 2.1*(10^9)])
xlim([0 101])
ylabel('Remaining energies in nJ')
xlabel('node number')
title(['Remaining energies after T3 cycles, T3=' num2str(T3) ' cycles - corner sink node'])