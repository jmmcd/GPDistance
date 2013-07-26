all:
	cd java && make
	cd python && make

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
	cd java && make sampleOneStepProbabilities
	cd python && make sampleOneStepProbabilities
