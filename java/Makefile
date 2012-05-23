cp=${CLASSPATH}:.:libraries/jgrapht-jdk1.5.jar:libraries/jscheme-7.2.jar

all:
	cd TreeDistance && make
	cd GP && make

mh:
	java -Xmx512m -cp ${cp} GP.Sample metro

completeMatrices:
	java -cp ${cp} GP.Sample completeMatrices

testdistances:
	java -cp ${cp} TreeDistance.OverlapDistance
	java -cp ${cp} TreeDistance.TreeDistance
