classdef N4BiasCorrection
    %N4BIASCORRECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = N4BiasCorrection(inputArg1,inputArg2)
            %N4BIASCORRECTION Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

% #include <SimpleITK.h>
% #include <iostream>
% #include <stdlib.h>
% 
% 
% namespace sitk = itk::simple;
% 
% int main ( int argc, char* argv[] ) {
% 
%   if ( argc < 2 ) {
%     std::cerr << "Usage: N4BiasFieldCorrection inputImage outputImage";
%     std::cerr << " [shrinkFactor] [maskImage] [numberOfIterations]";
%     std::cerr << " [numberOfFittingLevels]\n";
%     return 1;
%   }
% 
%   sitk::Image inputImage = sitk::ReadImage( argv[1], sitk::sitkFloat32 );
%   sitk::Image image = inputImage;
% 
%   sitk::Image maskImage;
%   if ( argc > 4 ) {
%     maskImage = sitk::ReadImage( argv[4], sitk::sitkUInt8 );
%   } else {
%     maskImage = sitk::OtsuThreshold( image, 0, 1, 200 );
%   }
% 
%   unsigned int shrinkFactor = 1;
%   if ( argc > 3 ) {
%     shrinkFactor = atoi( argv[3] );
%     std::vector<unsigned int> shrink( inputImage.GetDimension(), shrinkFactor );
%     image = sitk::Shrink( inputImage, shrink );
%     maskImage = sitk::Shrink( maskImage, shrink );
%   }
% 
%   sitk::N4BiasFieldCorrectionImageFilter *corrector
%     = new sitk::N4BiasFieldCorrectionImageFilter();
% 
%   unsigned int numFittingLevels = 4;
% 
%   if ( argc > 6) {
%     numFittingLevels = atoi(argv[6]);
%   }
% 
%   if ( argc > 5 ) {
%     unsigned int it = atoi( argv[5] );
%     std::vector<unsigned int> iterations( numFittingLevels, it );
%     corrector->SetMaximumNumberOfIterations( iterations );
%   }
% 
%   sitk::Image corrected_image = corrector->Execute( image, maskImage );
% 
%   sitk::Image log_bias_field = corrector->GetLogBiasFieldAsImage( inputImage );
% 
%   sitk::Image corrected_image_full_resolution = sitk::Divide( inputImage, sitk::Exp( log_bias_field ) );
% 
%   sitk::WriteImage( corrected_image_full_resolution, argv[2] );
% 
%   if (shrinkFactor > 1) {
%     sitk::WriteImage( corrected_image, "CXX-Example-N4BiasFieldCorrection-shrunk.nrrd" );
%   }
% 
%   if (getenv("SITK_NOSHOW") == NULL)
%     sitk::Show(corrected_image, "N4 Corrected");
% 
%   return 0;
% }