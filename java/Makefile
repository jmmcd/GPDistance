cp=${CLASSPATH}:.:libraries/jgrapht-jdk1.5.jar:libraries/jscheme-7.2.jar

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
	java -cp ${cp} GP.Sample uniformSampleMatrices 6
