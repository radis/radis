MatLab interface
---------------
RADIS can be accessed from Matlab by setting up Python environment in Matlab ::

    s = py.radis.calc_spectrum(1900,2300,molecule='CO',isotope='1,2,3',pressure=1.01325,Tgas=700,mole_fraction=0.1,path_length=1,databank='hitran');
    s.apply_slit(0.5,'nm');
    s.plot(show=true);
