function Cannnyedgedetector
close all;
%image reading and converting to Grayscale
image = 'image.png';
Original = imread(image);
figure(1);
subplot(1,2,1);
imshow(Original);
title("Original Image");
if(length(Original(1,1,:) )== 3)
     Original = Rgb2gs(Original); % Convert rgb to grayscale
end
subplot(1,2,2);
imshow(Original);
title('Original Grayscale Image');
%image reading and converting to Grayscale


% Gaussian Blur and Convolution
Kernel = Gauss(); % Calculate Kernel Matrix of Gaussian
Blurred = ApplyFilter(Original,Kernel); % Apply blurring
figure(2);
imshow(Blurred);
title('Blurred Image');
%Gaussian Blur and Convolution

KernelX = [1 0 -1; 1 0 -1; 1  0 -1]; % Prewitt Operator For X 
KernelY = [1 1  1; 0 0  0;-1 -1 -1]; % Prewitt Operator For Y

[DetectedMag,DetectedIntens,DetX,DetY] = getMagnitude(Blurred,KernelX,KernelY);

figure(3);
subplot(1,2,1);
imshow(DetX);
title("X axis");
subplot(1,2,2);
imshow(DetY);
title("Y axis");

figure(4);
imshow(DetectedMag);
title("Detected Diagonal");
% Calculating Magnitude and showing

EdgeThinned=EdgeThin(DetectedMag,DetectedIntens);


figure(5); 
imshow(EdgeThinned);
title("Edge Thinned");

Thresholded=Tresh(EdgeThinned);
figure(6);
imshow(Thresholded);
title("Thresholded");


end
function  Kernel = Gauss()
Kernel = zeros(5,5);
W = 0;
sigma = 1;
for i = 1:5
    for j = 1:5
        distance = (i - 3)^2 + (j - 3)^2;
        Kernel(i,j) = exp(-1 * (distance) / (2 * sigma * sigma));
        W = W + Kernel(i,j);
    end
end
Kernel = Kernel / W; % Normalization

end
function output = ApplyFilter(Matrix,Kernel)
[m,n] = size(Matrix); % size of image matrix
Matrix2 = Matrix; % copy image matrix
count = length(Kernel(1,:));
count = count-1; % count equals to kernel sub 1
for i =1 : m-count
    for j =1 :n-count
        temp = Matrix2(i:i+count , j: j+count); % Take 3 by 3 matrix
        temp =double(temp); % Convert to Double
        convolution = temp.*Kernel ; % multiply all elements with kernel Matrix
        Output(i,j) = sum(convolution(:)); % And sum multiplied matrix and save
    end
    output = uint8(Output);   
    
end
end
function ColorCode = Rgb2gs(Matrix)
    Red     = Matrix(:,:,1); % get Red color matrix of an image matrix    
    Green   = Matrix(:,:,2); % get Green color matrix of an image matrix    
    Blue    = Matrix(:,:,3); % get Blue color matrix of an image matrix    
    ColorCode = 0.2126 * Red + 0.7152 * Green + 0.0722 * Blue; % Calculation of grayscale
end
function [Detected,DetectedIntens, DetX,DetY]= getMagnitude(Matrix,KernelX,KernelY)
    [m,n] = size(Matrix);
    count =  length(KernelX(1,:));
    count = count-1;
    for i =1 : m-count
        for j =1 :n-count
            temp = double(Matrix(i:i+2,j:j+2)); % Save a 3 by 3 matrix to temporory value
            S1 = sum(sum(double(KernelX).* temp)); % Multiply values in a 3 by 3 matrix with Kernel X and Sum
            DetX(i,j) = uint8(S1); % convert 8 bit integer 
            S2 = sum(sum(double(KernelY).* temp)); % Multiply values in a 3 by 3 matrix with Kernel Y and Sum
            DetY(i,j) = uint8(S2); % convert 8 bit integer
            Detected(i,j) = sqrt( double(S1*S1 + S2*S2)); % Calculate magnitute like hypotenuse
            DetectedIntens(i,j) = atan2(S2,S1); % Calculate 2-dimensional arctanjant 
        end
    end
    Detected = Norm(Detected);
end
function Z = EdgeThin(GMatrix,IMatrix)
[m n] = size(GMatrix);
I = IMatrix*180./pi; % Convert from radian to degree
Z = uint8(zeros(m,n)); % create 8 bit integer zero matrix that size is m and n
for i=1:m
    for j=1:n
        if( I(i,j)< 0) 
            I(i,j)= 180 + I(i,j); % sum lower than 0 values with 180
        end
    end
end

for i=2: m-1
    for j=2: n-1
        q= 255;
        l= 255;
        if( 0 <= I(i,j)< 22.5 || 157.5 <= I(i,j)< 180) % to find pixel values that are neighborn by intensity
            q= GMatrix(i,j+1);
            l= GMatrix(i,j-1);
        elseif(22.5<= I(i,j)< 67.5)% to find pixel values that are neighborn by intensity
            q=GMatrix(i+1,j-1);
            l= GMatrix(i-1,j+1);
        elseif(67.5<= I(i,j)< 112.5)% to find pixel values that are neighborn by intensity
            q=GMatrix(i+1,j);
            l= GMatrix(i-1,j);
        elseif(112.5<= I(i,j)< 157.5)% to find pixel values that are neighborn by intensity
            q=GMatrix(i-1,j-1);
            l= GMatrix(i+1,j+1);
        elseif(GMatrix(i,j) >= q) && (GMatrix(i,j)>=l)% to find pixel values that are neighborn by intensity
                    Z(i,j) = GMatrix(i,j);
        end
        if(GMatrix(i,j)>=q && GMatrix(i,j)>=l) % Compare with intensity neighborns
            Z(i,j) = GMatrix(i,j);  % If neignborns have lower pixel, mark main pixel as an edge pixel
        else
            Z(i,j) =0;  % if both of them are  higher than main pixel, remove from the edge
        end
        
    end
    end
end
function [Matrix] = Tresh(Matrix)
HighRate = 0.09; % To calculate HighThreshold value
LowRate = 0.05;  % To calculate LowThreshold value
[m n] = size(Matrix);

HighThreshold = HighRate* max(max(double(Matrix))); % HighThreshold by multiplying max value of image and HighRate
LowThreshold =  HighThreshold * LowRate; % Multiply HighTreshold with LowRate
weak = uint8(25);
strong = uint8(255);
for i=1:m
    for j=1:n
        if(Matrix(i,j)>= HighThreshold)% If a pixel color value is higher than
            Matrix(i,j) = strong;      % HighThreshold, equal to strong value.
        elseif(Matrix(i,j) < LowThreshold)% If a pixel color value is lower than
            Matrix(i,j) = 0;              % HighThreshold, equal to 0 value.
        elseif(Matrix(i,j) < HighThreshold && Matrix(i,j)>LowThreshold) % If it is between 2 values
            Matrix(i,j) = weak;           % Equals to weak value
        end
    end
end
end
function Normalized = Norm(Matrix)
Max = max(max(Matrix));
Normalized= ( double(Matrix )* 255)/ Max;
Normalized = uint8(Normalized);
end
