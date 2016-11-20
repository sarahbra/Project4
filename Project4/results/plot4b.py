from matplotlib.pyplot import *

number_spins = 20.
f = open('results4c0T2.400000.txt')
f2 = open('results4c1T2.400000.txt')
mcs = []
mcs2 = []
E = []
E2 = []
Mabs = []
Mabs2 = []
acc_conf = []
acc_conf2 = []

for line in f.readlines()[1:]:
    temp = line.split(",")
    mcs.append(float(temp[0]))
    E.append(float(temp[1]))
    Mabs.append(float(temp[4]))
    acc_conf.append(float(temp[5]))
f.close()
for line in f2.readlines()[1:]:
    temp = line.split(",")
    mcs2.append(float(temp[0]))
    E2.append(float(temp[1]))
    Mabs2.append(float(temp[4]))
    acc_conf2.append(float(temp[5]))
f2.close()
E = np.array(E);
Mabs = np.array(Mabs);

E = E/number_spins**2
Mabs = Mabs/number_spins**2

E2 = np.array(E2);
Mabs2 = np.array(Mabs2);

E2 = E2/number_spins**2
Mabs2 = Mabs2/number_spins**2

steadystate = 20000
ordered = [5864,2128392]
unordered = [19969,2151366]
T = [1.0,2.4]

figure()
#subplot(2,1,1)

plot( T,unordered,label="Unordered")
plot(T,ordered,label="Ordered")
title("Accepted configurations of kbT")
#axis([0,10000])
xlabel("kbT")
ylabel("Accepted configurations(kbT)")
legend()

"""
subplot(2,1,2)
plot( mcs,Mabs,label="abs(M(t))")
#axis([0,10000])
xlabel("t(Monte Carlo cycles)")
ylabel("abs(M(t))")
legend()
"""
show()
    
