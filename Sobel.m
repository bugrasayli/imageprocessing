function main
close all;
image = 'image.png';  % Value of image
Original = imread(image); % Read the image and asign
figure(1);  
subplot(1,2,1);
imshow(Original);   % Show the image in figure
title("Original Image");
if(length(Original(1,1,:)))
   Original = rgb2grayscale(Original); % control if image color code
end

subplot(1,2,2);
imshow(Original);
title("Converted Grayscale if Yes");

% RGB or Grayscale
title("Original Image");
kernelX = [-1 0 +1;                    % Sobel Edge Detector X filter
           -2 0 +2 ; 
           -1 0 +1];
kernelY = [1 2 1;                   % Sobel Edge Detector Y filter
           0  0  0; 
          -1 -2 -1];


              
[Output1,Output2, DetectedMag] = ApplyFilter(Original,kernelX,kernelY); %apply filter and equal 2 values

figure(2);
subplot(1,2,1);
imshow(Output1);    % Show Edge detected image in X axis
title("X axis");

subplot(1,2,2);
imshow(Output2);
title("Y axis");    % Show Edge detected image in Y axis

figure(3);
imshow(DetectedMag);
title("Edge Detection");


a= Threshold (DetectedMag,70);
figure(4);
imshow(a);
title("Thresholded - Final Image");



figure(5);

subplot(2,2,1); 
imshow(Threshold(DetectedMag,25)); title("25 thresh Level");

subplot(2,2,2); 
imshow(Threshold(DetectedMag,50)); title("50 thresh Level");

subplot(2,2,3); 
imshow(Threshold(DetectedMag,100)); title("100 thresh Level");

subplot(2,2,4); 
imshow(Threshold(DetectedMag,150)); title("150 thresh Level");










end
function [Output1 Output2 DetectedMag] = ApplyFilter(Matrix,KernelX,KernelY)
[m,n] = size(Matrix); % Size of Matrix
Padded = zeros(m+2,n+2) ; % Create zero matrix 
Padded(2:m+1,2:n+1) = Matrix(:,:); % copy image matrix to zero matrix for padding
    for i =1 : m
        for j =1 :n
            temp = Padded(i:i+2 , j: j+2); % Take 3 by 3 matrix from image matrix
            temp =double(temp); % convert to double
            convolution1 = temp.*KernelX ; % multiply kernel with all values in matrix
            convolution2 = temp.*KernelY ; % multiply kernel with all values in matrix
            S1= sum(convolution1(:)); % sum all values in convoluted matrix
            S2=sum(convolution2(:));    % sum all values in convoluted matrix
            DetectedMag(i,j) = sqrt( (double(S1*S1)+ (double(S2*S2))));   % Calculate magnitude sqrt((X*X) + (Y*Y))
            Output1(i,j)=S1; % equal value to Output Matrix for X axis
            Output2(i,j)=S2;    % equal value to Output Matrix for Y axis
        end
    Output1 = uint8(Output1); %convert to uint(0 to 255)
    Output2 = uint8(Output2);    %convert to uint(0 to 255)
    end
    
    DetectedMag = uint8(DetectedMag);
    
end
function Matrix=rgb2grayscale(Matrix)
    R = Matrix(:,:,1); % Get Red color code of image matrix
    G = Matrix(:,:,1); % Get Green color code of image matrix
    B = Matrix(:,:,1); % Get Blue color code of image matrix
    Matrix  = (0.3 * R) + (0.59 * G) + (0.11 * B);
    % Weighted method to calculate grayscale 

end
function Thresholded = Threshold(Matrix,ThreshLevel)
    [m,n] = size(Matrix);
    for i=1:m
        for j=1:n
            if(Matrix(i,j) <= ThreshLevel) % control if pixels are lower than threshlevel
                Thresholded(i,j) = 0; %Remove
            else
                Thresholded(i,j) = Matrix(i,j); %Don’t do anything
            end
        end
    end
end
