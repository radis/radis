''' 
The AdBKurucz class below is inspired by an Exojax class developed by Hiroyuki Tako ISHIKAWA.
It allows loading data from the Kurucz database and performing several calculations on it.

Author: Racim Menasria
Date: June 2023
'''




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
from tqdm import tqdm
import os
from contextlib import closing
import pkgutil
import io
from io import BytesIO
import periodictable
 

class AdBKurucz():
    ccgs = 29979245800.0
    ecgs = 4.80320450e-10
    mecgs = 9.10938356e-28

    def __init__(self):
        self.kurucz_url_base = "http://kurucz.harvard.edu/linelists/gfall/gf"
        self.hdf5_file = None
        self.data = None
        self.pfTdat, self.pfdat = self.load_pf_Barklem2016()
        self.populations = None
        


    def get_url(self, atomic_number, ionization_state):
        ionization_state = str(ionization_state).zfill(2)
        return f"http://kurucz.harvard.edu/linelists/gfall/gf{atomic_number}{ionization_state}.all"
  

    def load_pf_Barklem2016(self):
        """load a table of the partition functions for 284 atomic species.

        * See Table 8 of Barklem & Collet (2016); https://doi.org/10.1051/0004-6361/201526961

        Returns:
            pfTdat (pd.DataFrame): steps of temperature (K)
            pfdat (pd.DataFrame): partition functions for 284 atomic species
        """
        pfT_str = 'T[K], 1.00000e-05, 1.00000e-04, 1.00000e-03, 1.00000e-02, 1.00000e-01, 1.50000e-01, 2.00000e-01, 3.00000e-01, 5.00000e-01, 7.00000e-01, 1.00000e+00, 1.30000e+00, 1.70000e+00, 2.00000e+00, 3.00000e+00, 5.00000e+00, 7.00000e+00, 1.00000e+01, 1.50000e+01, 2.00000e+01, 3.00000e+01, 5.00000e+01, 7.00000e+01, 1.00000e+02, 1.30000e+02, 1.70000e+02, 2.00000e+02, 2.50000e+02, 3.00000e+02, 5.00000e+02, 7.00000e+02, 1.00000e+03, 1.50000e+03, 2.00000e+03, 3.00000e+03, 4.00000e+03, 5.00000e+03, 6.00000e+03, 7.00000e+03, 8.00000e+03, 9.00000e+03, 1.00000e+04'
        pfTdat = pd.read_csv(io.StringIO(pfT_str), sep=',') 
        pfTdat = pd.Series(pfTdat.columns[1:]).astype('float64') # Converts the values to float64, skipping the first value

        # read the local file instead of the one in the exojax package
        with open('./pfdat.txt', 'r') as f:
            pfdat = pd.read_csv(f, sep='\s+', comment='#', names=pfTdat.index)

        #print(pfdat.head())
        #print("pfdat", pfdat)
        print("len(pfdat)", len(pfdat))
        print("len(pfTdat)", len(pfTdat))

        return pfTdat, pfdat



    def Sij0(self, A, gupper, nu_lines, elower, QTref_284, QTmask, Irwin=False):
        """Reference Line Strength in Tref=296K, S0.

    Args:
       A: Einstein coefficient (s-1)
       gupper: the upper state statistical weight
       nu_lines: line center wavenumber (cm-1)
       elower: elower
       QTref_284: partition function Q(Tref)
       QTmask: mask to identify a rows of QTref_284 to apply for each line
       Irwin: if True(1), the partition functions of Irwin1981 is used, otherwise those of Barklem&Collet2016

    Returns:
       Sij(T): Line strength (cm)
    """

    # Assign Q(Tref) for each line
        QTref = np.zeros_like(QTmask, dtype=float)
        for i, mask in enumerate(QTmask):
            QTref[i] = QTref_284[mask]

        # Use Irwin_1981 for Fe I (mask==76)  #test211013Tako
        if Irwin == True:
            QTref[np.where(QTmask == 76)[0]] = self.partfn_Fe(Tref_original)

        S0 = -A*gupper*np.exp(-hcperk*elower/Tref_original)*np.expm1(-hcperk*nu_lines/Tref_original)\
            / (8.0*np.pi*ccgs*nu_lines**2*QTref)

        return S0
    
    def load_atomicdata(self):
        """load atomic data and solar composition.

        * See  Asplund et al. 2009, Gerevesse et al. 1996

        Returns:
            ipccd (pd.DataFrame): table of atomic data

        Note:
            atomic.txt is in data/atom
        """
        ipccc = ('ielem', 'ionizationE1', 'dam1', 'dam2',
                'solarA', 'mass', 'ionizationE2')
        adata = pkgutil.get_data('exojax', 'data/atom/atomic.txt')
        ipccd = pd.read_csv(BytesIO(adata), sep='\s+', skiprows=1,
                            usecols=[1, 2, 3, 4, 5, 6, 7], names=ipccc)
        return ipccd


    def load_ionization_energies(self):
        """Load atomic ionization energies.

        Returns:
            df_ionE (pd.DataFrame): table of ionization energies

        Note:
            NIST_Atomic_Ionization_Energies.txt is in data/atom
        """
        fn_IonE = pkgutil.get_data(
            'exojax', 'data/atom/NIST_Atomic_Ionization_Energies.txt')
        df_ionE = pd.read_csv(BytesIO(fn_IonE), sep='|', skiprows=6, header=0)
        return df_ionE

    def pick_ionE(self,ielem, iion, df_ionE):
        """Pick up ionization energy of a specific atomic species.

        Args:
            ielem (int): atomic number (e.g., Fe=26)
            iion (int): ionized level (e.g., neutral=1, singly ionized=2, etc.)
            df_ionE (pd.DataFrame): table of ionization energies

        Returns:
            ionE (float): ionization energy

        Note:
            NIST_Atomic_Ionization_Energies.txt is in data/atom
        """
        def f_droppare(x): return x.str.replace('(', '', regex=True).str.replace(')', '', regex=True).str.replace(
            '[', '', regex=True).str.replace(']', '', regex=True).str.replace('                                      ', '0', regex=True)
        ionE = float(f_droppare(df_ionE[(df_ionE['At. num '] == ielem) & (
            df_ionE[' Ion Charge '] == iion-1)]['      Ionization Energy (a) (eV)      ']))
        return ionE
    
    def download_file(self):
        """Download a file from an url to a specified output path."""
        
        # extracts the file's name from l'URL
        filename = self.url.split("/")[-1]
        
        with closing(requests.get(self.url, stream=True)) as r:
            total_size = int(r.headers.get('content-length', 0))
            block_size = 1024
            with open(filename, 'wb') as f, tqdm(
                    total=total_size, 
                    unit='B', 
                    unit_scale=True, 
                    unit_divisor=1024, 
                    desc=filename
                ) as pbar:
                for data in r.iter_content(block_size):
                    f.write(data)
                    pbar.update(len(data))
        return filename



    def air_to_vac(self,wlair):
        """Convert wavelengths [AA] in air into those in vacuum.

        * See http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

        Args:
            wlair:  wavelengthe in air [Angstrom]
            n:  Refractive Index in dry air at 1 atm pressure and 15ºC with 0.045% CO2 by volume (Birch and Downs, 1994, Metrologia, 31, 315)

        Returns:
            wlvac:  wavelength in vacuum [Angstrom]
        """
        s = 1e4 / wlair
        n = 1. + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s*s) + \
            0.0001599740894897 / (38.92568793293 - s*s)
        wlvac = wlair * n
        return wlvac

    def read_kurucz(self,kuruczf):
        """Input Kurucz line list (http://kurucz.harvard.edu/linelists/)

        Args:
            kuruczf: file path

        Returns:
            A:  Einstein coefficient in [s-1]
            nu_lines:  transition waveNUMBER in [cm-1] (#NOT frequency in [s-1])
            elower: lower excitation potential [cm-1] (#converted from eV)
            eupper: upper excitation potential [cm-1] (#converted from eV)
            gupper: upper statistical weight
            jlower: lower J (rotational quantum number, total angular momentum)
            jupper: upper J
            ielem:  atomic number (e.g., Fe=26)
            iion:  ionized level (e.g., neutral=1, singly)
            gamRad: log of gamma of radiation damping (s-1) #(https://www.astro.uu.se/valdwiki/Vald3Format)
            gamSta: log of gamma of Stark damping (s-1)
            gamvdW:  log of (van der Waals damping constant / neutral hydrogen number) (s-1)
        """
        with open(kuruczf) as f:
            lines = f.readlines()
        wlnmair, loggf, species, elower, jlower, labellower, eupper, jupper, labelupper, \
            gamRad, gamSta, gamvdW, ref, \
            NLTElower, NLTEupper, isonum, hyperfrac, isonumdi, isofrac, \
            hypershiftlower, hypershiftupper, hyperFlower, hypernotelower, hyperFupper, hypternoteupper, \
            strenclass, auto, landeglower, landegupper, isoshiftmA \
            = \
            np.zeros(len(lines)), np.zeros(len(lines)), np.array(['']*len(lines), dtype=object), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), \
            np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), \
            np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), \
            np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines)), \
            np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(
                len(lines)), np.zeros(len(lines)), np.zeros(len(lines))
        ielem, iion = np.zeros(len(lines), dtype=int), np.zeros(
            len(lines), dtype=int)

        for i, line in enumerate(lines):
            wlnmair[i] = float(line[0:11])
            loggf[i] = float(line[11:18])
            species[i] = str(line[18:24])
            ielem[i] = int(species[i].split('.')[0])
            iion[i] = int(species[i].split('.')[1])+1
            elower[i] = float(line[24:36])
            jlower[i] = float(line[36:41])
            eupper[i] = float(line[52:64])
            jupper[i] = float(line[64:69])
            gamRad[i] = float(line[80:86])
            gamSta[i] = float(line[86:92])
            gamvdW[i] = float(line[92:98])

        elower_inverted = np.where((eupper-elower) > 0,  elower,  eupper)
        eupper_inverted = np.where((eupper-elower) > 0,  eupper,  elower)
        jlower_inverted = np.where((eupper-elower) > 0,  jlower,  jupper)
        jupper_inverted = np.where((eupper-elower) > 0,  jupper,  jlower)
        elower = elower_inverted
        eupper = eupper_inverted
        jlower = jlower_inverted
        jupper = jupper_inverted

        wlaa = np.where(wlnmair < 200, wlnmair*10, self.air_to_vac(wlnmair*10))
        nu_lines = 1e8 / wlaa[::-1]  # [cm-1]<-[AA]
        loggf = loggf[::-1]
        ielem = ielem[::-1]
        iion = iion[::-1]
        elower = elower[::-1]
        eupper = eupper[::-1]
        jlower = jlower[::-1]
        jupper = jupper[::-1]
        glower = jlower*2+1
        gupper = jupper*2+1
        A = 10**loggf / gupper * (self.ccgs*nu_lines)**2 \
            * (8*np.pi**2*self.ecgs**2) / (self.mecgs*self.ccgs**3)
        gamRad = gamRad[::-1]
        gamSta = gamSta[::-1]
        gamvdW = gamvdW[::-1]


        data_dict = {
                "A": A,
                "nu_lines": nu_lines,
                "elower": elower,
                "eupper": eupper,
                "gupper": gupper,
                "jlower": jlower,
                "jupper": jupper,
                "ielem": ielem,
                "iion": iion,
                "gamRad": gamRad,
                "gamSta": gamSta,
                "gamvdW": gamvdW,
            }

        self.data = pd.DataFrame(data_dict)
        return self.data



        

    def store_hdf5(self,data, output_path):
        """Store data in a HDF5 file."""

        #df = pd.DataFrame(data)
        data.to_hdf(output_path, key='kurucz_data')

    def read_hdf5(self,file_path):
        """Read data from a HDF5 file."""

        return pd.read_hdf(file_path)

    def partfn(self, key, T):
        """Partition function from Barklem & Collet (2016).

        Args:
        atom: the atom name
        T: temperature

        Returns:
        partition function Q
        """
        print(f"Température: {T}")
        # Locate the row for the specific atom and ionization state
        pf_atom = self.pfdat.loc[f'{key}']
        # Extract temperature and partition function values
        pfT_values = self.pfTdat.values.flatten()[1:]  # Exclude the first value (it's the temperature unit)
        pf_values = pf_atom.values[1:]  # Exclude the first value (it's the atomic number)

        pfT_values = pfT_values.astype(float)
        pf_values = pf_values.astype(float)

        #print("pfT_values:", pfT_values)
        #print("pf_values:", pf_values)


        # Interpolate to find the partition function at the desired temperature
        Q = np.interp(T, pfT_values, pf_values)
        #print(f"Partition function: {Q}")

        return Q

        
    def calculate_populations(self, atom, temperature):
        # Select the partition function for the specific atom
        #atom_pf = self.pfdat.loc[f'{atom}_I']  

        # Calculate the partition function at a certain temperature
        QT_atom = self.partfn(atom, temperature)

        # Calculate energy/temperature ratio
        energy_temp_ratio = self.data["elower"] / (0.695 * temperature)
        #print(f'Energy/temp ratio: {energy_temp_ratio}')   # print the energy to temperature ratio for debugging

        # Calculate level populations using Boltzmann statistics
        self.populations = np.exp(-energy_temp_ratio) / QT_atom

        return self.populations




    def plot_spectrum(self,data, populations, temperature):
        intensities = data['A'] * populations 
        #print(f"Intensities: {intensities}")  
        plt.figure(figsize=(10, 6))
        plt.plot(data['nu_lines'], intensities)
        plt.title(f'Spectrum at {temperature}K')
        plt.xlabel('Wave Number (cm-1)')
        plt.ylabel('Intensity')
        plt.show()

    def process (self, atomic_number, ionization_state, temperature):
        self.kurucz_file = f"gf{atomic_number}{ionization_state}.all"
        self.hdf5_file = f"gf{atomic_number}{ionization_state}.hdf5"
        self.url = self.get_url(atomic_number, ionization_state) 
        self.kuruczf = self.download_file()  
        self.data = self.read_kurucz(self.kuruczf)
        self.store_hdf5(self.data,self.hdf5_file)
        self.data = self.read_hdf5(self.hdf5_file)
        
        # Convert atomic_number to element symbol
        element_symbol = periodictable.elements[int(atomic_number)].symbol

        # Construct the key

        if ionization_state == '00':
            key = element_symbol + '_I'
        elif ionization_state =='01':
            key = element_symbol + '_II'
        else :
            key= element_symbol + '_III'

        populations=self.calculate_populations(key,temperature)
        self.plot_spectrum(self.data, populations, temperature)