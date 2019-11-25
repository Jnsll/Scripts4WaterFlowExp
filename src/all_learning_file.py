import LearningFile as lf

folder = "/run/media/jnsll/b0417344-c572-4bf5-ac10-c2021d205749/exps_modflops/results/"
indicator = "H"

for site in range(1,9):
    print("site", site)
    lf.createCSVFileForASite(site, folder)

for site in range(10,31):
    print("site", site)
    lf.createCSVFileForASite(site, folder)