clc;
clear all;
close all;


%Constant Preparation
HEIGHT = 400;
WIDTH = 400;
CONST_A = 0.001;
CONST_B = 1;
CONST_P = 0.5;
LINE_NUM = 26;


% read image A and B
A = imread('1.jpg');
B = imread('2.jpg');
%read image A
imshow(A);
hold on;
for i = 1:LINE_NUM
    [xA(i,1), yA(i, 1)] = ginput(1);
    plot(xA(i,1),yA(i,1),'.','color',[1 0 0]); 
    [xA(i,2), yA(i, 2)] = ginput(1);
    plot(xA(i,2),yA(i,2),'.','color',[1 0 0]); 
    plot(xA(i,:), yA(i,:),'-rs');
end
hold off;
%read image B
imshow(B);
hold on;
for i = 1:LINE_NUM
    [xB(i,1), yB(i, 1)] = ginput(1);
    plot(xB(i,1),yB(i,1),'.','color',[1 0 0]); 
    [xB(i,2), yB(i, 2)] = ginput(1);
    plot(xB(i,2),yB(i,2),'.','color',[1 0 0]); 
    plot(xB(i,:), yB(i,:),'-rs');
end
hold off;


cons = 1;
for medicons = 0:0.01:1
% median of every pair of lines of A and B
xM = (1-medicons)*xA + medicons*xB;
yM = (1-medicons)*yA + medicons*yB;
%Perpendicular Matrix
perpendM = [0 1; -1 0];


%%%%%%AAAAAAAAAAAAAA%%%%%%
% line PQ
p = [xM(:,1) yM(:,1)];
q = [xM(:,2) yM(:,2)];
lines = q - p;
% line PQ prime
pp = [xA(:,1) yA(:,1)];
qp = [xA(:,2) yA(:,2)];
linesp = qp - pp;
%distort A
distortedA(1:HEIGHT,1:WIDTH,1:3) = 255;
%dis AND weight
for x = 1:WIDTH
    for y = 1:HEIGHT
        pix = [x y];
        xp = [x-p(:,1) y-p(:,2)];
        perpendlines=lines*perpendM;
        normlines = dot(lines, lines, 2).^0.5;
        u = dot(xp, lines, 2)./(normlines.^2);
        v = dot(xp, perpendlines, 2)./normlines;
        perpendlinesp = linesp * perpendM;
        normlinesp = dot(linesp, linesp, 2).^0.5;
        pixp = pp + [linesp(:,1).*u  linesp(:,2).*u] + [perpendlinesp(:,1) .* v ./normlinesp perpendlinesp(:,2) .* v ./normlinesp];
        d = [pixp(:,1)-pix(1) pixp(:,2)-pix(2)];
        dSum = [0 0];
        weightSum = 0;
        for i = 1:LINE_NUM
            % dist
            if u(i)<0
                dist = dot(p(i,:)-pix,p(i,:)-pix,2).^0.5;
            elseif u(i) > 1
                dist = dot(q(i,:)-pix,q(i,:)-pix,2).^0.5;
            else
                dist = abs(v(i));
            end
            weight = ((dot(lines(i,:),lines(i,:),2)^0.5)^CONST_P/(CONST_A + dist))^CONST_B;
            weightSum = weightSum + weight; 
            dSum = dSum + d(i,:) * weight;
        end
        pixPrime = int64(pix + dSum/weightSum);
        if(pixPrime(1) > 0 && pixPrime(1) <= WIDTH && pixPrime(2) > 0 && pixPrime(2) <= HEIGHT)
            distortedA(y, x, :) = A(pixPrime(2),pixPrime(1),:);
        end
    end
end


%%%%%%BBBBBBBBBBBBBB%%%%%%
% line PQ
p = [xM(:,1) yM(:,1)];
q = [xM(:,2) yM(:,2)];
lines = q - p;
% line PQ prime
pp = [xB(:,1) yB(:,1)];
qp = [xB(:,2) yB(:,2)];
linesp = qp - pp;
%distort B
distortedB(1:HEIGHT,1:WIDTH,1:3) = 255;
%dis AND weight
for x = 1:WIDTH
    for y = 1:HEIGHT
        pix = [x y];
        xp = [x-p(:,1) y-p(:,2)];
        perpendlines=lines*perpendM;
        normlines = dot(lines, lines, 2).^0.5;
        u = dot(xp, lines, 2)./(normlines.^2);
        v = dot(xp, perpendlines, 2)./normlines;
        perpendlinesp = linesp * perpendM;
        normlinesp = dot(linesp, linesp, 2).^0.5;
        pixp = pp + [linesp(:,1).*u  linesp(:,2).*u] + [perpendlinesp(:,1) .* v ./normlinesp perpendlinesp(:,2) .* v ./normlinesp];
        d = [pixp(:,1)-pix(1) pixp(:,2)-pix(2)];
        dSum = [0 0];
        weightSum = 0;
        for i = 1:LINE_NUM
            % dist
            if u(i)<0
                dist = dot(p(i,:), pix, 2).^0.5;
            elseif u(i) > 1
                dist = dot(q(i,:), pix, 2).^0.5;
            else
                dist = abs(v(i));
            end
            weight = ((dot(lines(i,:),lines(i,:),2)^0.5)^CONST_P/(CONST_A + dist))^CONST_B;
            weightSum = weightSum + weight; 
            dSum = dSum + d(i,:) * weight;
        end
        pixPrime = int64(pix + dSum/weightSum);
        if(pixPrime(1) > 0 && pixPrime(1) <= WIDTH && pixPrime(2) > 0 && pixPrime(2) <= HEIGHT)
            distortedB(y, x, :) = B(pixPrime(2),pixPrime(1),:);
        end
    end
end


%show the morphing
imwrite(uint8((1-medicons)*distortedA+medicons*distortedB),strcat('new',num2str(cons),'.jpg'));
cons = cons+1;
end