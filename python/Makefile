all:
	echo "Nothing to do here."

completeMatricesDepth1:
	cd RandomWalks && ./random_walks.py ../../results/depth_1/

completeMatricesDepth2:
	cd RandomWalks && ./random_walks.py ../../results/depth_2/

completeMatricesGALength4::
	cd RandomWalks && ./random_walks.py ../../results/ga_length_4/

compareTPCalculatedVSampledDepth1:
	cd RandomWalks && ./plotting.py compareTPCalculatedVSampled ../../results/depth_1

compareTPCalculatedVSampledDepth2:
	cd RandomWalks && ./plotting.py compareTPCalculatedVSampled ../../results/depth_2

compareTPCalculatedVSampledGALength4:
	cd RandomWalks && ./plotting.py compareTPCalculatedVSampled ../../results/ga_length_4

writeSteadyStateDepth1:
	cd RandomWalks && ./plotting.py writeSteadyState ../../results/depth_1

writeSteadyStateDepth2:
	cd RandomWalks && ./plotting.py writeSteadyState ../../results/depth_2

writeSteadyStateGALength4:
	cd RandomWalks && ./plotting.py writeSteadyState ../../results/ga_length_4

makeGridPlotsDepth1:
	cd RandomWalks && ./plotting.py makeGridPlots ../../results/depth_1

makeGridPlotsDepth2:
	cd RandomWalks && ./plotting.py makeGridPlots ../../results/depth_2

makeGridPlotsGALength4:
	cd RandomWalks && ./plotting.py makeGridPlots ../../results/ga_length_4

makeCorrelationTableDepth1:
	cd RandomWalks && ./plotting.py makeCorrelationTable ../../results/depth_1 "GP Depth 1"

makeCorrelationTableDepth2:
	cd RandomWalks && ./plotting.py makeCorrelationTable ../../results/depth_2 "GP Depth 2"

makeCorrelationTableGALength4:
	cd RandomWalks && ./plotting.py makeCorrelationTable ../../results/ga_length_4 "GA Length 4"

uniformSampleDepth6:
	cd RandomWalks && ./random_walks.py graphDistances ../../results/uniform_depth_6/

rwSampleDepth6:
	cd RandomWalks && ./random_walks.py graphDistances ../../results/rw_sample_depth_6/

mhSampleDepth6:
	cd RandomWalks && ./random_walks.py graphDistances ../../results/mh_sample_depth_6/

