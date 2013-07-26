cp=${CLASSPATH}:.:libraries/jgrapht-jdk1.6.jar:libraries/jscheme-7.2.jar:libraries/RTED_v1.1.jar

all:
	cd TreeDistance && make
	cd GP && make

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
	java -cp ${cp} GP.AllTrees 1

completeSpaceDepth2:
	java -cp ${cp} GP.AllTrees 2

uniformSampleDepth6:
	java -cp ${cp} GP.Sample uniformSampleMatrices 6 1000

rwSampleDepth6:
	java -cp ${cp} GP.Sample rwSampleMatrices 6

mhSampleDepth6:
	java -cp ${cp} GP.Sample mhSampleMatrices 6

randomWalkingDepth6:
	java -cp ${cp} GP.Sample randomWalking 6

sampleOneStepProbabilities:
	java -cp ${cp} GP.Sample sampleOneStepProbabilities 1
