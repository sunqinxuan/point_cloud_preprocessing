% clear all;close all;clc;

% ps=load('43_pointCloud.txt');
% im=im2double(imread('43_saveImage.jpg'));

for row=1:30
    for col=1:40
        gridpoint=zeros(256,7);
        first=row*40*1200+col*16;
        num=1;
        for m=1:16
            for n=1:16
                k=first+640*m+n;
                if ps(k,1)==0&&ps(k,2)==0&&ps(k,3)==0
                    continue;
                else
                    gridpoint(num,1)=0;
                    gridpoint(num,2)=ps(k,1);
                    gridpoint(num,3)=ps(k,2);
                    gridpoint(num,4)=ps(k,3);
                    gridpoint(num,5)=0;
                    gridpoint(num,6)=0;
                    gridpoint(num,7)=0;
                end
            end
        end
    end
end