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
sys.path.append("")
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

# Grafico fit e studio quantitativo bontà adattamento
fig, ax = plt.subplots(2,1, figsize=(9,6), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
fig.subplots_adjust(hspace=0)
ax[0].set_title('Fit istogramma considerando scattering di Rayleigh')
ax[0].hist(lam,bins=n_bins, weights=n_photons,color="limegreen", label="Dati" )
ax[0].plot(bincenters, y_fit, color="blue", label="Fit")
ax[0].set_ylabel('Numero fotoni')
ax[0].legend(frameon=False)
scarto=n[mask]-y_fit[mask]
ax[1].scatter(bincenters[mask],scarto, marker=".",color="darkred")
ax[1].set_xlabel("Lunghezza d'onda (nm)")
ax[1].set_ylabel('Scarti')      
plt.show()

#Chiqudro
chi2_2 =  np.sum( (scarto)**2 /n[mask] ) 
#Chi quadro ridotto
ndof_2 = len(n[mask])-len(par)
chi2_rid_2=chi2_2/ndof_2
print("Il chi quadro ridotto è: {:.3f}".format(chi2_rid_2))

print("La temperatura stimata è di {:.1f} +- {:.1f} K" .format(par[1],mt.sqrt(pcov[1][1])))