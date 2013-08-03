all:
	cd java && make
	cd python && make

depth1:
	mkdir -p results/depth_1
	cd java && make completeSpaceDepth1
	cd java && make completeMatricesDepth1
	cd python && make completeMatricesDepth1
	cd java && make sampleOneStepProbabilitiesDepth1
	cd python && make writeSteadyStateDepth1
	cd python && make makeCorrelationTableDepth1
	cd python && make makeGridPlotsDepth1
	cd python && make compareTPCalculatedVSampledDepth1

depth2:
	mkdir -p results/depth_2
	cd java && make completeSpaceDepth2
	cd java && make completeMatricesDepth2
	cd python && make completeMatricesDepth2
	cd java && make sampleOneStepProbabilitiesDepth2
	cd python && make writeSteadyStateDepth2
	cd python && make makeCorrelationTableDepth2
	cd python && make makeGridPlotsDepth2
	cd python && make compareTPCalculatedVSampledDepth2

ga_length4:
	mkdir -p results/ga_length_4
	cd python && make completeMatricesGALength4
	cd python && make writeSteadyStateGALength4
	cd python && make makeCorrelationTableGALength4
	cd python && make makeGridPlotsGALength4

completeMatricesDepth1:
	cd java && make completeMatricesDepth1
	cd python && make completeMatricesDepth1

completeMatricesDepth2:
	cd java && make completeMatricesDepth2
	cd python && make completeMatricesDepth2

uniformSampleMatricesDepth6:
	cd java && make uniformSampleDepth6
	cd python && make uniformSampleDepth6

rwSampleMatricesDepth6:
	cd java && make rwSampleDepth6
	cd python && make rwSampleDepth6

mhSampleMatricesDepth6:
	cd java && make mhSampleDepth6
	cd python && make mhSampleDepth6

sampleOneStepProbabilities:
	cd java && make sampleOneStepProbabilitiesDepth1
