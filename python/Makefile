all:
	echo "Nothing to do here."

completeMatricesDepth1:
	cd RandomWalks && ./random_walks.py ../../results/depth_1/

completeMatricesDepth2:
	cd RandomWalks && ./random_walks.py ../../results/depth_2/

compareTPCalculatedVSampledDepth1:
	cd RandomWalks && ./plotting.py compareTPCalculatedVSampled ../../results/depth_1 > ../../results/depth_1/compare_TP_calculated_v_sampled.txt

compareTPCalculatedVSampledDepth2:
	cd RandomWalks && ./plotting.py compareTPCalculatedVSampled ../../results/depth_2 > ../../results/depth_2/compare_TP_calculated_v_sampled.txt

writeSteadyStateDepth1:
	cd RandomWalks && ./plotting.py writeSteadyState ../../results/depth_1

writeSteadyStateDepth2:
	cd RandomWalks && ./plotting.py writeSteadyState ../../results/depth_2

makeGridPlotsDepth1:
	cd RandomWalks && ./plotting.py makeGridPlots ../../results/depth_1

makeGridPlotsDepth2:
	cd RandomWalks && ./plotting.py makeGridPlots ../../results/depth_2

makeCorrelationTableDepth1:
	cd RandomWalks && ./plotting.py makeCorrelationTable ../../results/depth_1 "GP Depth 1"

makeCorrelationTableDepth2:
	cd RandomWalks && ./plotting.py makeCorrelationTable ../../results/depth_2 "GP Depth 2"

uniformSampleDepth6:
	cd RandomWalks && ./random_walks.py ../../results/uniform_depth_6/

rwSampleDepth6:
	cd RandomWalks && ./random_walks.py ../../results/rw_sample_depth_6/

mhSampleDepth6:
	cd RandomWalks && ./random_walks.py ../../results/mh_sample_depth_6/

