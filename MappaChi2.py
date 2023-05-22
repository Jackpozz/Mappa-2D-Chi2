import numpy as np
import csv
import array as arr
from scipy import optimize
import matplotlib.pyplot as plt


with open("Risonanze.csv",'r') as i:                #Apre un file .csv con dati disposti come x,y,ErrY diviso da virgole
    rawdata = list(csv.reader(i,delimiter=","))   
    
exampledata = np.array(rawdata[0:],dtype=float)     #Crea 3 arrays contenenti i dati in input
xdata = exampledata[:,0]
ydata = exampledata[:,1]
sydata = exampledata[:,2]





def f(x,A,Ω,Δ):                                             #Definisce una funzione fantasma da usare quando si variano
    return A*x/np.sqrt(x**4-2*(Ω**2-2*Δ**2)*x**2+Ω**4)      #le coppie tra A,Ω,Δ (Non fa niente in realtà SIUM) 

num=0                                                       #Dichiara variabili e arrays riempiti di 0 da sostituire
ChiMin=0
pDeltaMin=0
pOmegaMin=0
Delta=arr.array('f',[])
Omega=arr.array('f',[])
Chi=[[0 for col in range(100)] for row in range(100)]       #Il valore in range(N) definisce la lunghezza di righe e colonne
AMin=[[0 for col in range(100)] for row in range(100)]      #creando una matrice virtuale NxN 
                                                            #(non aumentare troppo N o esplode il pc)


for i in range(8900,9500,6):                              #Ciclo per riempire due arrays di valori di prova (Toys)
    Delta.append(i)                                         #Le variabili in range() determinano rispettivamente la zona di
for j in range(481000,482000,10):                           #scansione Min, Max e lo Step.
    Omega.append(j)                                         #Di conseguenza (Max-Min)/Step deve restituire esattamente N precedente

   
for k in range(len(Omega)):                                                             #Per ogni coppia di Omega e Delta trova il
    for l in range(len(Delta)):                                                         #BestFit per A e crea una matrice di Chi2
        def g(x,A):
            return A*x/np.sqrt(x**4-2*(Omega[l]**2-2*Delta[k]**2)*x**2+Omega[l]**4)
        p_init = [6000] 
        p_best, pcov, infodict, mesg, ier = optimize.curve_fit(g, xdata, ydata, 
        sigma = sydata, p0=p_init, bounds=(0, +np.inf), full_output=True)
        
        AMin[k][l] = p_best[0]
        
        for m in range(len(xdata)):
            Chi[k][l] += (((g(xdata[m], p_best[0])-ydata[m]) / sydata[m])**2)
            
            
        if ((Chi[k][l] <= ChiMin) or ((l == 0) and (k == 0))):
            ChiMin = Chi[k][l]
            pDeltaMin = k
            pOmegaMin = l
            

ErrChi = Chi.copy()                                                         #Trova i valori di Omega e Delta per cui ΔChi<1
ErrX = arr.array('i',[])                                                    
ErrY = arr.array('i',[])
for n in range(len(Omega)):
    if (ErrChi[pOmegaMin][n] - ChiMin <= 1):
       ErrX.append(n)
       
    if (ErrChi[n][pDeltaMin] - ChiMin <= 1):
       ErrY.append(n)
       
    
       
b = ErrX[0]
c = ErrY[0]
d = ErrX[-1]
e = ErrY[-1]
ErrOmega = Omega[pOmegaMin]-Omega[b]
ErrDelta = Delta[pDeltaMin]-Delta[c]
ErrA = np.sqrt(pcov[0])   
            

print("Delta_Min (Hz) = ")
print(Delta[pDeltaMin], " +- ", ErrDelta)
print("Omega_Min (Hz) = ")
print(Omega[pOmegaMin], " +- ", ErrOmega)
print("A_Min = ")
print(AMin[pDeltaMin][pOmegaMin], " +- ", ErrA)
print("Chi_Min=")
print(Chi[pDeltaMin][pOmegaMin])
print("Matrice del Chi Quadro (NxN): ")
print(np.matrix(Chi))


plt.figure(1, figsize=(8,4.5))
plt.subplots_adjust(left=0.09, bottom=0.09, top=0.97, right=0.99)

image = plt.imshow(Chi, origin = 'lower', vmax = 11, cmap ='jet',                 #Plotta il grafico inserendo la matrice Chi
                   extent = [Omega[0], Omega[99], Delta[0], Delta[99]],           #vmax come valore di Chi massimo (si cambia per
                   aspect = (Omega[99]-Omega[0])/(Delta[99]-Delta[0]))            #aumentare il gradiente del colore)

clb = plt.colorbar(image) 
clb.ax.set_title("Chi Quadro") 
plt.xlim(Omega[0], Omega[99])
plt.ylim(Delta[0], Delta[99])
plt.xlabel(r'$Omega (kHz)$', fontsize=20)
plt.ylabel(r'$Delta (kHz)$', fontsize=20)
plt.hlines(Delta[e], Omega[0], Omega[99], linestyles = "dotted")
plt.hlines(Delta[c], Omega[0], Omega[99], linestyles = "dotted")
plt.vlines(Omega[b], Delta[0], Delta[99], linestyles = "dotted")
plt.vlines(Omega[d], Delta[0], Delta[99], linestyles = "dotted")
plt.savefig('Chivs3Var.png')
plt.show()
