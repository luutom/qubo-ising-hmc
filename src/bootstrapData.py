import sys
import numpy as np

data = open(sys.argv[1],"r").readlines()

m=[]
e=[]
si=[[] for i in range(7)]


for dat in data:
    if "#" not in dat:
        observables = dat.split()
        m.append(float(observables[3]))
        e.append(float(observables[4]))
        for j in range(7):
            si[j].append(float(observables[j+5]))


eMean = np.mean(np.array(e))
mMean = np.mean(np.array(m))
sMean=[]
for i in range(7):
    sMean.append(np.mean(np.array(si[i])))

NBS = 100
BSLENGTH=len(m)

sError=[]
bootstrapIndices = np.random.randint(0, BSLENGTH, [NBS, BSLENGTH])
mError = np.std(np.array([np.mean(np.array([m[cfg] for cfg in bootstrapIndices[sample]]), axis=0) for sample in range(NBS) ] ), axis=0)
eError = np.std(np.array([np.mean(np.array([e[cfg] for cfg in bootstrapIndices[sample]]), axis=0) for sample in range(NBS) ] ), axis=0)
for i in range(7):
    sError.append(np.std(np.array([np.mean(np.array([si[i][cfg] for cfg in bootstrapIndices[sample]]), axis=0) for sample in range(NBS) ] ), axis=0))

print("# <E>     <dE>       <m>        <dm>")
print(eMean," ",eError," ",mMean," ",mError)
print("")
print("# i   <S_i>     <dS_i>")
for i in range(7):
    print(i+1," ",sMean[i]," ",sError[i])
                  
