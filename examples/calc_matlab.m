//Example of using RADIS from matlab
//All methods of RADIS can be accessed using py.radis followed by the method to be accessed


s = py.radis.calc_spectrum(1900,2300,molecule='CO',isotope='1,2,3',pressure=1.01325,Tgas=700,mole_fraction=0.1,path_length=1,databank='hitran');
s.apply_slit(0.5,'nm');
s.plot(show=true);
