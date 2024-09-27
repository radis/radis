from radis import calc_spectrum
s = calc_spectrum(1900, 2300,         # cm-1
                  molecule='CO',
                  isotope='1,2,3',
                  pressure=1.01325,   # bar
                  Tgas=700,           # K
                  mole_fraction=0.1,
                  path_length=1,      # cm
                  databank='hitran',  # or 'hitemp', 'geisa', 'exomol'
                  )
s.apply_slit(0.5, 'nm')       # simulate an experimental slit
s.plot('radiance')