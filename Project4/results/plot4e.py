from matplotlib.pyplot import *
f40 = open('results4e_Spins40_mcs300000.txt')
T40 = []
E40 = []
Mabs40 = []
C_v40 = []
chi40 = []

for line in f40.readlines()[1:]:
    temp = line.split(",")
    T40.append(float(temp[2]))
    E40.append(float(temp[3]))
    Mabs40.append(float(temp[4]))
    C_v40.append(float(temp[5]))
    chi40.append(float(temp[6]))

f40.close()

f60 = open('results4e_Spins60_mcs500000.txt')
T60 = []
E60 = []
Mabs60 = []
C_v60 = []
chi60 = []

for line in f60.readlines()[1:]:
    temp = line.split(",")
    T60.append(float(temp[2]))
    E60.append(float(temp[3]))
    Mabs60.append(float(temp[4]))
    C_v60.append(float(temp[5]))
    chi60.append(float(temp[6]))

f60.close()

f100 = open('results4e_Spins100_mcs500000.txt')
T100 = []
E100 = []
Mabs100 = []
C_v100 = []
chi100 = []

for line in f100.readlines()[1:]:
    temp = line.split(",")
    T100.append(float(temp[2]))
    E100.append(float(temp[3]))
    Mabs100.append(float(temp[4]))
    C_v100.append(float(temp[5]))
    chi100.append(float(temp[6]))

f100.close()

f140 = open('results4e_Spins140_mcs700000.txt')
T140 = []
E140 = []
Mabs140 = []
C_v140 = []
chi140 = []

for line in f140.readlines()[1:]:
    temp = line.split(",")
    T140.append(float(temp[2]))
    E140.append(float(temp[3]))
    Mabs140.append(float(temp[4]))
    C_v140.append(float(temp[5]))
    chi140.append(float(temp[6]))

f140.close()

figure()

plot(T40,E40,label="40*40")
plot(T60,E60,label="60*60")
plot(T100,E100,label="100*100")
plot(T140,E140, label="140*140")
title("Energy with respect to kbT")
xlabel("kbT")
ylabel("E(kbT)")
legend()
show()

figure()

plot(T40,Mabs40,label="40*40")
plot(T60,Mabs60,label="60*60")
plot(T100,Mabs100,label="100*100")
plot(T140,Mabs140, label="140*140")
title("Absolute value of M with respect to kbT")
xlabel("kbT")
ylabel("abs(M(kbT))")
legend()
show()

figure()

plot(T40,C_v40,label="40*40")
plot(T60,C_v60,label="60*60")
plot(T100,C_v100,label="100*100")
plot(T140,C_v140, label="140*140")
title("Heat capacity with respect to kbT")
xlabel("kbT")
ylabel("C_v(kbT)")
#legend()
show()

figure()

plot(T40,chi40,label="40*40")
plot(T60,chi60,label="60*60")
plot(T100,chi100,label="100*100")
plot(T140,chi140, label="140*140")
title("Magnetic susceptibility with respect to kbT")
xlabel("kbT")
ylabel("chi(kbT)")
legend()
show()
