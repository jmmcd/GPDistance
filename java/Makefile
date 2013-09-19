cp=${CLASSPATH}:.:libraries/jgrapht-jdk1.6.jar:libraries/jscheme-7.2.jar:libraries/RTED_v1.1.jar

all:
	cd TreeDistance && make
	cd GP && make

testTP:
	java -Xmx512m -cp ${cp} GP.Mutation "(+ x (+ y y))" "y"

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

mhSampleDepth6:
	java -cp ${cp} GP.Sample mhSampleMatrices 6

randomWalkingDepth2:
	java -cp ${cp} GP.Sample randomWalkingMatrices 2

sampleForSuperNodeDepth2:
	java -cp ${cp} GP.Sample sampleForSuperNode 2

randomWalkingDepth6:
	java -cp ${cp} GP.Sample randomWalkingMatrices 6

sampleOneStepProbabilitiesDepth1:
	java -cp ${cp} GP.Sample sampleOneStepProbabilities 1

sampleOneStepProbabilitiesDepth2:
	java -cp ${cp} GP.Sample sampleOneStepProbabilities 2
