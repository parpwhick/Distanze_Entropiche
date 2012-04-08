function make()

	fprintf('Compiling mexfiles ...');
	
    mex -largeArrayDims -lgomp computeAbsPower.cpp -outdir ../clustering/
    mex -largeArrayDims -lgomp getSparseDerivativeMatrix.cpp -outdir ../clustering/
    mex -largeArrayDims -lgomp pNormPowParallel.cpp -outdir ../clustering/
    mex -largeArrayDims -lgomp pNormPowDegParallel.cpp -outdir ../clustering/

	fprintf('\nThe mexfiles have been moved to the folder ../clustering/\n');

end