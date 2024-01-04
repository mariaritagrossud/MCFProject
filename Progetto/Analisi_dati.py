"""
File con analisi dati sulla stella X 
_________________________________________
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd
import math as mt
from scipy.optimize import curve_fit
import sys

#Importo modulo con funzioni
sys.path.append("/Users/mariarita/Metodi_computazionali/gitHub/MCF/Progetto ")
import func as f

#Carico file con dati
data = pd.read_csv("observed_starX.csv", sep=",")
lam = data["lambda (nm)"].array
n_photons = data["photons"].array

#Visualizzo distribuzione
n_bins = 100
n,bins, p = plt.hist(lam,bins=n_bins, weights=n_photons, color="limegreen") 
plt.title("Istogramma distribuzione fotoni in funzione lunghezza d'onda")
plt.ylabel("Numero fotoni")
plt.xlabel("Lunghezza d'onda (nm)")
plt.show()

#Faccio fit e plotto risultati

bincenters = (bins[:-1] + bins[1:])/2
mask = np.nonzero(n)
par, pcov = curve_fit(f.n_photons_fit, xdata=bincenters[mask], ydata=n[mask], sigma=np.sqrt(n[mask]), 
                    p0=[3.309509356824346e-30,8001], absolute_sigma=True)

y_fit = f.n_photons_fit(bincenters,par[0],par[1])

plt.figure(figsize=(7,6))
n,bins, p = plt.hist(lam,bins=n_bins, weights=n_photons,color="limegreen", label="Dati") 
plt.plot(bincenters, y_fit, color="blue", label="Fit")
plt.title("Istogramma distribuzione fotoni in funzione lunghezza d'onda")
plt.ylabel("Numero fotoni")
plt.xlabel("Lunghezza d'onda (nm)")
plt.legend(frameon=False)
plt.show()
print("La temperatura stimata è di {:.1f} +- {:.1f} K" .format(par[1],mt.sqrt(pcov[1][1])))