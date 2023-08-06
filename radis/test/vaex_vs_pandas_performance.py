import time

import matplotlib.pyplot as plt

from radis import calc_spectrum


# Compares the time performance of pandas and Vaex
def compare_vaex_pandas_time():
    time_list, timeC_list, lines_list = [], [], []
    time_list_va, timeC_list_va, lines_list_va = [], [], []
    wmin = 1000

    for engine in ["vaex", "pandas"]:
        for w_range in [
            10,
            30,
            100,
            150,
            500,
            600,
            700,
            800,
            900,
            1000,
            1500,
            3000,
        ]:  # [100, 200, 500, 900]:
            t0 = time.time()
            s = calc_spectrum(
                wmin,
                wmin + w_range,  # cm-1
                molecule="H2O",
                isotope="1,2,3",
                pressure=1.01325,  # bar
                Tgas=1000,
                mole_fraction=0.1,
                databank="hitemp",  # or 'hitemp'
                diluent="air",
                verbose=0,
                engine=engine,
            )
            t1 = time.time()
            if engine == "vaex":
                timeC_list_va.append(s.conditions["calculation_time"])
                lines_list_va.append(s.conditions["lines_calculated"])
                time_list_va.append(t1 - t0)
                # lines_list_va.append(s.conditions['lines_calculated']+s.conditions['lines_cutoff'])
            else:
                timeC_list.append(s.conditions["calculation_time"])
                lines_list.append(s.conditions["lines_calculated"])
                time_list.append(t1 - t0)
                # lines_list.append(s.conditions['lines_calculated']+s.conditions['lines_cutoff'])

    plt.figure()
    plt.plot(lines_list, time_list, "k", label="pandas total")
    plt.plot(lines_list, timeC_list, "k--", label="pandas computation")
    plt.plot(lines_list_va, time_list_va, "r", label="vaex total")
    plt.plot(lines_list_va, timeC_list_va, "r--", label="vaex computation")
    plt.ylabel("Time [s]")
    plt.xlabel("Number of lines")
    plt.legend()


# Compare the memory performance of Pandas and Vaex
def compare_pandas_vs_vaex_memory():

    import tracemalloc

    # Change range of wav to obtain result for different number of lines , repeat same procedure both vaex and pandas
    # Loop is not used because it gives better results for memory performance for vaex or pandas which one is used second time
    tracemalloc.start()
    s = calc_spectrum(
        1000,
        2800,  # cm-1
        molecule="H2O",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=1000,  # K
        mole_fraction=0.1,
        wstep="auto",
        path_length=1,  # cm
        databank="hitemp",  # or 'hitemp', 'geisa', 'exomol'
        engine="vaex",
    )

    s.apply_slit(0.5, "nm")  # simulate an experimental slit
    memory = tracemalloc.get_traced_memory()
    print((memory[1]))
    tracemalloc.stop()
    print(s.conditions["lines_calculated"])
