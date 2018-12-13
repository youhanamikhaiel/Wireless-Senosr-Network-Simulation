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
kbits = Data + Overhead;
d0 = sqrt(eps_short/eps_long);
xy = randi([0 n],2,n);
sink = [w/2; l/2];
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l])
xlim([0 w])
hold on
scatter(sink(1),sink(2),'filled')
title('Network topology')

%Plotting the range of nodes with single and dual hops
R = 30;
angle = 0:0.01:2*pi;
xc = R*cos(angle)+sink(1);
yc = R*sin(angle)+sink(2);
hold on
plot(xc,yc)

%Calculating the distance between each node and the sink
dist = [];
for i=1:n
    dist(i) = sqrt(((xy(1,i)-sink(1))^2) + ((xy(2,i)-sink(2))^2));
end

% Initializing arrays to define the nodes' death
Ei1 = Ei;
death = [];
logic_dead = ones(1,n);
ii = [];
dead_sign = ones(1,n);

% Starting the transmission of the nodes 
for j=1:4000
    for i=1:n
        m = 0;
        energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2) + kbits*Eagg;
        energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4) + kbits*Eagg;
        energy_loss_rec = kbits*Eelec;
        
        if dist(i) <= R
            if dist(i) <= d0
                if energy_loss_short < Ei1(i)
                    Ei1(i) = Ei1(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            else
                if energy_loss_long < Ei1(i)
                    Ei1(i) = Ei1(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            end
        else
            [distmin1, distmin2, ii(j,i)] = mindistance(xy,i,sink,logic_dead,R);
            
            energy_loss_short2 = kbits*Eelec + kbits*eps_short*(distmin1^2) + kbits*Eagg;
            energy_loss_short3 = kbits*Eelec + kbits*eps_short*(distmin2^2) + kbits*Eagg;
            if ii(j,i) ~= 0
                if energy_loss_short2 < Ei1(i)
                    if energy_loss_short3 < Ei1(ii(j,i))
                        Ei1(i) = Ei1(i) - energy_loss_short2;
                        Ei1(ii(j,i)) = Ei1(ii(j,i)) - energy_loss_short3;
                    else
                        logic_dead(ii(j,i)) = 0;
                        if dead_sign(i) == 1;
                            death(ii(j,i)) = j;
                            dead_sign(ii(j,i)) = 0;
                    end
                    end
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            end
        end
    end
end

T1 = min(death)
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
hold on 
plot(ones(n+1)*T1,0:n)

%%
Ei2 = Ei;
death = [];
logic_dead = ones(1,n);
ii = [];
dead_sign = ones(1,n);

for j=1:T1
    for i=1:n
        m = 0;
        energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2) + kbits*Eagg;
        energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4) + kbits*Eagg;
        energy_loss_rec = kbits*Eelec;
        
        if dist(i) <= R
            if dist(i) <= d0
                if energy_loss_short < Ei2(i)
                    Ei2(i) = Ei2(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            else
                if energy_loss_long < Ei2(i)
                    Ei2(i) = Ei2(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            end
        else
            [distmin1, distmin2, ii(j,i)] = mindistance(xy,i,sink,logic_dead,R);
            
            energy_loss_short2 = kbits*Eelec + kbits*eps_short*(distmin1^2) + kbits*Eagg;
            energy_loss_short3 = kbits*Eelec + kbits*eps_short*(distmin2^2) + kbits*Eagg;
            if ii(j,i) ~= 0
                if energy_loss_short2 < Ei2(i)
                    if energy_loss_short3 < Ei2(ii(j,i))
                        Ei2(i) = Ei2(i) - energy_loss_short2;
                        Ei2(ii(j,i)) = Ei2(ii(j,i)) - energy_loss_short3;
                    else
                        logic_dead(ii(j,i)) = 0;
                        if dead_sign(i) == 1;
                            death(ii(j,i)) = j;
                            dead_sign(ii(j,i)) = 0;
                    end
                    end
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            end
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


%%
Ei3 = Ei;
death = [];
logic_dead = ones(1,n);
ii = [];
dead_sign = ones(1,n);

% Starting the transmission of the nodes
for RR=10:10:80
    for j=1:4000
        for i=1:n
            m = 0;
            energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2) + kbits*Eagg;
            energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4) + kbits*Eagg;
            energy_loss_rec = kbits*Eelec;
            
            if dist(i) <= RR
                if dist(i) <= d0
                    if energy_loss_short < Ei3(i)
                        Ei3(i) = Ei3(i) - energy_loss_short;
                    else
                        logic_dead(i) = 0;
                        if dead_sign(i) == 1;
                            death(i) = j;
                            dead_sign(i) = 0;
                        end
                    end
                else
                    if energy_loss_long < Ei3(i)
                        Ei3(i) = Ei3(i) - energy_loss_short;
                    else
                        logic_dead(i) = 0;
                        if dead_sign(i) == 1;
                            death(i) = j;
                            dead_sign(i) = 0;
                        end
                    end
                end
            else
                [distmin1, distmin2, ii(j,i)] = mindistance(xy,i,sink,logic_dead,RR);
                
                energy_loss_short2 = kbits*Eelec + kbits*eps_short*(distmin1^2) + kbits*Eagg;
                energy_loss_short3 = kbits*Eelec + kbits*eps_short*(distmin2^2) + kbits*Eagg;
                if ii(j,i) ~= 0
                    if energy_loss_short2 < Ei3(i)
                        if energy_loss_short3 < Ei3(ii(j,i))
                            Ei3(i) = Ei3(i) - energy_loss_short2;
                            Ei3(ii(j,i)) = Ei3(ii(j,i)) - energy_loss_short3;
                        else
                            logic_dead(ii(j,i)) = 0;
                            if dead_sign(i) == 1;
                                death(ii(j,i)) = j;
                                dead_sign(ii(j,i)) = 0;
                            end
                        end
                    else
                        logic_dead(i) = 0;
                        if dead_sign(i) == 1;
                            death(i) = j;
                            dead_sign(i) = 0;
                        end
                    end
                end
            end
        end
    end
    T11(RR/10) = min(death);
end
RR = 10:10:80;
figure
plot (RR,T11)
xlim([0 90])
ylim([0 max(T11)+2])
xlabel('Maximum distance for direct transmission')
ylabel('Number of cycles to the death of the first node T1')
title('T1 versus R')
Ropt = RR(find(T11==max(T11)));

%%
Ei4 = Ei;
death = [];
logic_dead = ones(1,n);
ii = [];
dead_sign = ones(1,n);

for j=1:max(T11)
    for i=1:n
        m = 0;
        energy_loss_short = kbits*Eelec + kbits*eps_short*(dist(i)^2) + kbits*Eagg;
        energy_loss_long = kbits*Eelec + kbits*eps_long*(dist(i)^4) + kbits*Eagg;
        energy_loss_rec = kbits*Eelec;
        
        if dist(i) <= Ropt(end)
            if dist(i) <= d0
                if energy_loss_short < Ei4(i)
                    Ei4(i) = Ei4(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            else
                if energy_loss_long < Ei4(i)
                    Ei4(i) = Ei4(i) - energy_loss_short;
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
            end
        else
            [distmin1, distmin2, ii(j,i)] = mindistance(xy,i,sink,logic_dead,R);
            
            energy_loss_short2 = kbits*Eelec + kbits*eps_short*(distmin1^2) + kbits*Eagg;
            energy_loss_short3 = kbits*Eelec + kbits*eps_short*(distmin2^2) + kbits*Eagg;
            if ii(j,i) ~= 0
                if energy_loss_short2 < Ei4(i)
                    if energy_loss_short3 < Ei4(ii(j,i))
                        Ei4(i) = Ei4(i) - energy_loss_short2;
                        Ei4(ii(j,i)) = Ei4(ii(j,i)) - energy_loss_short3;
                    else
                        logic_dead(ii(j,i)) = 0;
                        if dead_sign(i) == 1;
                            death(ii(j,i)) = j;
                            dead_sign(ii(j,i)) = 0;
                    end
                    end
                else
                    logic_dead(i) = 0;
                    if dead_sign(i) == 1;
                        death(i) = j;
                        dead_sign(i) = 0;
                    end
                end
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
title('Remaining energies after T1 cycle for optimum R')

%%
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
kbits = Data + Overhead;
d0 = sqrt(eps_short/eps_long);
xy = randi([0 n],2,n);
sink = [w/2; l+125];
figure
scatter(xy(1,:),xy(2,:),'filled')
ylim([0 l+250])
xlim([-125 w+125])
hold on
scatter(sink(1),sink(2),'filled')
title('Network topology with far sink')
dist22 = [];
for i=1:n
    dist22(i) = sqrt(((xy(1,i)-sink(1))^2) + ((xy(2,i)-sink(2))^2));
end

farnodes = find(dist22<=132);
scatter(xy(1,farnodes),xy(2,farnodes),'filled')

%Plotting the range of nodes with single and dual hops
for R=10:10:130;
    angle = 0:0.01:2*pi;
    xc = R*cos(angle)+sink(1);
    yc = R*sin(angle)+sink(2);
    hold on
    plot(xc,yc)
end

stem(Ei)
title('Remaining energies for far sink with R < 125')
