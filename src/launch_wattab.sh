for RATE in 2.0 7.0 30.0 90.0 182.0 365.0 730.0 3652.0
do
	python3 watertab.py -site 11 -chr 0 -approx 0 -rate $RATE -f /run/media/jnsll/b0417344-c572-4bf5-ac10-c2021d205749/exps_modflops/results
done
