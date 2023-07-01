from radis.api.kuruczapi import AdBKurucz  

def main():
    # Create an instance of the AdBKurucz class
    spectrum = AdBKurucz()

    # Select element and temperature
    # for instance 2600 stands for Fe and 1101 for Na+
    atomic_number = '26'
    ionization_state = '00'
    temperature = 5000

    # Process with the data and plot the spectrum
    spectrum.process(atomic_number, ionization_state, temperature)

if __name__ == "__main__":
    main()