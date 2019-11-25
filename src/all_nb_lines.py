import GlobalIndicatorFileCreation as gifc

folder = "/run/media/jnsll/b0417344-c572-4bf5-ac10-c2021d205749/exps_modflops/results/"
indicator = "H"

for site in range(1,8):
    print("site", site)
    for chronicle in range(2):
        for approx in range(2):
            gifc.createCSVFileForASiteAndAContext(indicator, folder, site, chronicle, approx)


for site in range(10,31):
    print("site", site)
    for chronicle in range(2):
        for approx in range(2):
            gifc.createCSVFileForASiteAndAContext(indicator, folder, site, chronicle, approx)

