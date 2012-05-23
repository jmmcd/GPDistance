all:
	cd java && make
	cd python && make

completeMatrices:
	cd java && make completeMatrices
	cd python && make completeMatrices
