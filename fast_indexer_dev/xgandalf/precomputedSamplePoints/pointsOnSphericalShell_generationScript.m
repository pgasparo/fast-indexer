% needs "S2 Sampling Toolbox" from http://de.mathworks.com/matlabcentral/fileexchange/37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere/all_files
% needs "Rotation Matrix 3D" from https://de.mathworks.com/matlabcentral/fileexchange/46419-rotation-matrix-3d
thisFilesPath = fileparts(mfilename('fullpath'));

pitches = [0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
tolerances = [0, 0.005, 0.01, 0.02, 0.03, 0.04];

rotationMatrix1 = RotationMatrix([pi/8, pi/10, pi/12], 'eulerAngles').GetRotationMatrix();
rotationMatrix2 = RotationMatrix([pi/12, pi/8, pi/10], 'eulerAngles').GetRotationMatrix();
A = 4*pi;
for pitch = pitches
    pointCountOnSphere = round(A/pitch^2);
    for tolerance = tolerances

        symmetricAuxiliarySpheresCount = ceil(tolerance/pitch);
        radialPitch = tolerance/symmetricAuxiliarySpheresCount;
        if(pointCountOnSphere < 1300)   %very expensive operation
            samplePoints = ParticleSampleSphere('N', pointCountOnSphere, 'Nitr', 1.3e3); 
        else
            samplePoints = SpiralSampleSphere(pointCountOnSphere);
        end
        
        for i = 1:symmetricAuxiliarySpheresCount
            adaptedPointCount = pointCountOnSphere*(1-0.5*i/symmetricAuxiliarySpheresCount);
            if(pointCountOnSphere < 1300)  %very expensive operation
                tmp = ParticleSampleSphere('N', adaptedPointCount, 'Nitr', 1.3e3); 
            else
                tmp = SpiralSampleSphere(adaptedPointCount);
            end
            samplePoints = [samplePoints ; (rotationMatrix1^i*tmp')'*(1-i*radialPitch) ; (rotationMatrix2^i*tmp')'*(1+i*radialPitch)];
        end
        samplePoints(samplePoints(:,3) < 0,:) = [];


        filename = [thisFilesPath , '/pitch', num2str(pitch), '_tolerance', num2str(tolerance)];
        dlmwrite(filename, samplePoints,'delimiter','\t');

    end
end

filename = [thisFilesPath , '/pitches'];
dlmwrite(filename, pitches,'delimiter','\t');
filename = [thisFilesPath , '/tolerances'];
dlmwrite(filename, tolerances,'delimiter','\t');