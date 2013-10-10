cp=${CLASSPATH}:.:libraries/jscheme-7.2.jar:libraries/RTED_v1.1.jar

all:
	cd TreeDistance && make
	cd GP && make

testTP:
	java -Xmx512m -cp ${cp} GP.Mutation "(+ y (* x x))" "y"

testRandomWalk:
	java -Xmx512m -cp ${cp} GP.Sample testRandomWalk 4

mh:
	java -Xmx512m -cp ${cp} GP.Sample metro

completeMatricesDepth1:
	java -cp ${cp} GP.Sample completeMatrices 1

completeMatricesDepth2:
	java -cp ${cp} GP.Sample completeMatrices 2

testdistances:
	java -cp ${cp} TreeDistance.OverlapDistance
	java -cp ${cp} TreeDistance.TreeDistance

completeSpaceDepth1:
	java -cp ${cp} GP.AllTrees 1 | grep -v java | grep -v Number > ../results/depth_1/all_trees.dat

completeSpaceDepth2:
	java -cp ${cp} GP.AllTrees 2 | grep -v java | grep -v Number > ../results/depth_2/all_trees.dat

uniformSampleDepth6:
	java -cp ${cp} GP.Sample uniformSampleMatrices 6 1000

rwSampleDepth6:
	java -cp ${cp} GP.Sample rwSampleMatrices 6

rwSampleDepth2:
	java -cp ${cp} GP.Sample rwSampleMatrices 2

mhSampleDepth6:
	java -cp ${cp} GP.Sample mhSampleMatrices 6

randomWalkingDepth1:
	java -cp ${cp} GP.Sample randomWalking 1

randomWalkingDepth2:
	java -cp ${cp} GP.Sample randomWalking 2

randomWalkingDepth6:
	java -cp ${cp} GP.Sample randomWalking 6

sampleForSuperNodeDepth2:
	java -cp ${cp} GP.Sample sampleForSuperNode 2

sampleForSuperNodeDepth6:
	java -cp ${cp} GP.Sample sampleForSuperNode 6

sampleOneStepProbabilitiesDepth1:
	java -cp ${cp} GP.Sample sampleOneStepProbabilities 1

sampleOneStepProbabilitiesDepth2:
	java -cp ${cp} GP.Sample sampleOneStepProbabilities 2
