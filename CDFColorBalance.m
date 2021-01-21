function Main
close all;
filename = 'image.jpg';
global Red; global Green; global Blue;
[Red Green Blue] = GetRGB(filename);
[Hue Saturation Value] =Convert2HSV(Red,Green,Blue);
% 
% [Hue2 Saturation2 Value2] = rgb2hsv (Red,Green,Blue);

% [Hue2 Satur2 Val2] = rgb2hsv(Red,Green,Blue); 
% For testing conversation rgb to hsv with Matlab Function

%Calculate histogram values
[HueHist,SaturationHist,ValueHist]=calcHisofHSV(Hue,Saturation,Value);

%Draw Histogram graphs
DrawHistogramWithBar(HueHist,SaturationHist,ValueHist,'Original Hue','Original Saturation','Original Value');

%Calculate CDF Values
[HueHist,SaturationHist, ValueHist] = CalculateCumulative(HueHist,SaturationHist,ValueHist);

%Draw CDF
DrawCDF(HueHist,SaturationHist,ValueHist,'Original Hue CDF','Original Saturation CDF','Original Value CDF');

%Calculate Equalization
[HueHist,SaturationHist,ValueHist] = CalculateEqualization(HueHist,SaturationHist,ValueHist);

%Equalite Values
[Hue Saturation Value]=NewHSV(Hue,Saturation,Value,HueHist,SaturationHist,ValueHist);

%Calculate Last Histogram Values
[HueHist_New SaturationHist_New ValueHist_New] = calcHisofHSV(Hue,Saturation,Value);

%Draw last Histogram
DrawHistogramWithBar(HueHist_New,SaturationHist_New,ValueHist_New,'Last Hue Hist','Last Saturation Hist','Last Value Hist');

%Calculate Cumulative Values
[HueHist_New, ValueHist_New,SaturationHist_New] = CalculateCumulative(HueHist_New,SaturationHist_New,ValueHist_New);

%Draw CDF
DrawCDF(HueHist_New,SaturationHist_New,ValueHist_New,'Hue Last CDF','Saturation Last CDF','Value Last CDF');

[New_red New_green New_blue] = GetRGBfromHSV(Hue,Saturation,Value);

end
function [Red Green Blue] = GetRGB(image)
im = imread(image);
figure;
imshow(im);
title('Image');
Red = im(:,:,1);
Green = im(:,:,2);
Blue = im(:,:,3);
end
function [HueHist SaturationHist ValueHist]=calcHisofHSV(Hue,Saturation,Value)
%Color Matrix for Hue, Saturation and Value
N = length(Hue(:,1));
M = length(Hue(1,:));
HueHist = zeros(100,5);
SaturationHist = zeros(100,5);
ValueHist = zeros(100,5);

%Create color arrays that store color values from 0 to 100 
for i = 1: 101
   HueHist(i,1) = i-1;
   SaturationHist(i,1) = i-1;
   ValueHist(i,1) = i-1;
end

for i = 1: N
    for j = 1: M
        % Hue Color is temporary value for Image(i,j) color code
        % First index of array is 1, so sum with 1.
        Hue_ColorCode = floor(Hue(i,j));
        HueHist(Hue_ColorCode+1,2) = HueHist(Hue_ColorCode+1,2) + 1; 
        
        Sat_Colorcode = floor(Saturation(i,j)); 
        SaturationHist(Sat_Colorcode+1,2) = SaturationHist(Sat_Colorcode+1,2) +1; 
        
        Value_Colorcode = floor(Value(i,j)); 
        ValueHist(Value_Colorcode+1,2) = ValueHist(Value_Colorcode+1,2) +1; 
    end
end
end
function [Hue,Saturation,Value] = CalculateCumulative(Hue,Saturation,Value)
        global Red;
        N = length(Red(1,:));
        M = length(Red(:,1));
        %Calculate first color values.
        Hue(1,3)=Hue(1,2)/(N*M); 
        Saturation(1,3)=Saturation(1,2)/(N*M);
        Value(1,3)=Value(1,2) / (N * M);
        
        Hue(1,4)        = Hue(1,2) ;
        Saturation(1,4) = Saturation(1,2);
        Value(1,4)      = Value(1,2) ;
        
for i=2 : 101
    
        %CDF values / Sum with last value.
        Hue(i,3)        =Hue(i,2) / (N*M) + Hue(i-1,3); 
        Saturation(i,3) = Saturation(i,2)/(N*M) + Saturation(i-1,3);
        Value(i,3)      = Value(i,2) / (N * M) + Value(i-1,3);
        
        %Cumulative Histogram values
         Hue(i,4)        = Hue(i,2) + Hue(i-1,4); 
         Saturation(i,4) = Saturation(i,2) + Saturation(i-1,4);
         Value(i,4)      = Value(i,2) + Value(i-1,4);
end
end
function DrawHistogramWithBar(hue,saturation,value,HueName,SaturationName,ValueName)
figure();
set(gcf, 'Name', 'Histogram', 'NumberTitle', 'Off','position', [500, 0, 400, 800]);

subplot(3,1,1);
bar([hue(:,1) hue(:,2)],3,'FaceColor','y','EdgeColor',[0 0 0],'LineWidth',1.5);
title(HueName);

subplot(3,1,2);
bar([saturation(:,1) saturation(:,2)],3,'FaceColor','y','EdgeColor',[0 0 0],'LineWidth',1.5);
title(SaturationName);

subplot(3,1,3);
bar([value(:,1) value(:,2)],3,'FaceColor','y','EdgeColor',[0 0 0],'LineWidth',1.5);
title(ValueName);

end
function DrawCDF(Hue,Saturation,Value,HueName,SaturationName,ValueName)
figure();
set(gcf, 'Name', 'CDF', 'NumberTitle', 'Off','position', [500, 0, 400, 800]);

subplot(3,1,1);
plot(Hue(:,3))
title(HueName)

subplot(3,1,2);
plot(Saturation(:,3 ));
title(SaturationName);

subplot(3,1,3);
plot(Value(:,3));
title(ValueName);

end
function [red green blue] = GetRGBfromHSV(Hue,Saturation,Value)
Hue = Hue/100;
Saturation=Saturation/100;
Value = Value / 100;

HSV(:,:,1) =Hue;
HSV(:,:,2) =Saturation;
HSV(:,:,3) = Value;

figure();
rgb = hsv2rgb(HSV);

red = rgb(:,:,1);
green = rgb(:,:,2);
blue = rgb(:,:,3);
imshow(rgb);
end
function [HueHist ValueHist SaturationHist] = CalculateEqualization (HueHist,ValueHist,SaturationHist)
% 0 to 100 means 101 number of color code
ColorNumber =101;
L_Value = ColorNumber - 1;
for i=1:101
HueHist(i,5) = floor(L_Value * HueHist(i,3));
ValueHist(i,5) = floor(L_Value * ValueHist(i,3));
SaturationHist(i,5) = floor(L_Value * SaturationHist(i,3));
end
end
function [Hue,Saturation,Value]=NewHSV(Hue,Saturation,Value,HueHist,SaturationHist,ValueHist)
global Red;
N = length(Red(:,1));
M = length(Red(1,:));
% Changing last color values with new values that came 
% from Histogram Equalization
for i=1:N
    for j=1:M
        Hue(i,j)       = HueHist( floor(Hue(i,j) + 1),5); 
        Saturation(i,j)= SaturationHist( floor((Saturation(i,j)) + 1) ,5); 
        Value(i,j)     = ValueHist(floor((Value(i,j)) + 1) ,5); 
    end
end
end
function [hue sat val]=Convert2HSV(red, green, blue)
[hue sat val] = rgb2hsv(red,green,blue);
hue = hue * 100;
sat = sat * 100;
val = val * 100;

end
  

