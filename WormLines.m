%%
close all

%% Turn Worm to XY cordinates
%Worm is black (0) and background is white (1)
%Assume length 1 between each pixel
%buttom left corner is [0,0]
%x direction is horizontal left->right, y direction is vertical down-up
ymax=size(Worm,2);
xmax=size(Worm,1);
[rows,cols]=find(~Worm);
x=flipud(cols);
y=flipud(rows);

%Test
Fig1=figure;
Ax1=axes(Fig1);
hold(Ax1,'on');
imshow(Worm,'parent',Ax1);
scatter(x,y,0.5,'r','filled','parent',Ax1);

%% find prime direction of worm t and starting points
P=[x,y];
m=numel(x); %amount of points
coeff=pca(P);
t=coeff(:,1); %prime direction in column vec
n=coeff(:,2); %normal direction in column vec

CG=mean(P); %row vec
Q=zeros(m,2); %project points initalized
for k=1:m
    Q(k,:)=ProjectPointsOnLine(P(k,:),CG,t');
end
[~,LeftEdgeInd]=min((Q-CG)*t);
[~,RightEdgeInd]=max((Q-CG)*t);
E_l=P(LeftEdgeInd,:);
E_r=P(RightEdgeInd,:);
%projecting all points P on line t and taking the
%most extreme one

%Test
Fig2=figure;
Ax2=axes(Fig2);
hold(Ax2,'on');
imshow(Worm,'parent',Ax2);
scatter(Ax2,E_l(1),E_l(2),20,'b','filled');
scatter(Ax2,E_r(1),E_r(2),20,'b','filled');
scatter(Ax2,CG(1),CG(2),20,'r','marker','*');
quiver(Ax2,E_l(1),E_l(2),E_r(1)-E_l(1),E_r(2)-E_l(2));

%% Split worm to N segments in the prime direction
%use a different technique, rotation first.
N=30; %arbitrary
P_tn=(coeff*(P-CG)')'; %coeff is rotation matrix so t will become the x axis.
[~,LeftEdgeInd_tn]=min(P_tn(:,1));
[~,RightEdgeInd_tn]=max(P_tn(:,1));
E_l_tn=P_tn(LeftEdgeInd_tn,:);
E_r_tn=P_tn(RightEdgeInd_tn,:);

z=linspace(E_l_tn(1),E_r_tn(1),N);
P_SecMean_tn=zeros(N-1,2);
% FigTemp=figure;
% AxTemp=axes(FigTemp);
% hold(AxTemp,'on');
for k=2:N
   P_SecInd_tn=find(P_tn(:,1)>z(k-1) & P_tn(:,1)<z(k));
   P_Sec_tn=P_tn(P_SecInd_tn,:);
   P_SecMean_tn(k-1,:)=mean(P_Sec_tn,1);
%    scatter(AxTemp,P_Sec_tn(:,1),P_Sec_tn(:,2),10,'b','marker','*');
%    scatter(AxTemp,P_SecMean_tn(k-1,1),P_SecMean_tn(k-1,2),20,'r','marker','d');
%    pause(0.1);
end

%Test
Fig3=figure;
Ax3=axes(Fig3);
hold(Ax3,'on');
scatter(P_tn(:,1),P_tn(:,2),2,'r','filled','parent',Ax3);
scatter(Ax3,P_SecMean_tn(:,1),P_SecMean_tn(:,2),20,'marker','*');
% scatter(Ax3,E_l_tn(1),E_l_tn(2),20,'b','filled');
% scatter(Ax3,E_r_tn(1),E_r_tn(2),20,'b','filled');
for k=1:N
    zlines_x=[z(k),z(k)];
    plot(Ax3,zlines_x,Ax3.YLim,'linewidth',0.2,'linestyle','--','color',[0.5,0.5,0.5]);
end
%% calculate curvature from the segmented points on tn plane
Npeaks=3;
Px_tn=(smooth(P_SecMean_tn(:,1)));
Py_tn=(smooth(P_SecMean_tn(:,2)));
Curvature=NumCurvature(Px_tn,Py_tn);
Curvature_RGB=Curvature2RGB(Curvature);

Curv_Threshold=mean(Curvature);
[pks,locs]=findpeaks(Curvature,'MinPeakHeight',Curv_Threshold,...
    'SortStr','ascend');
Pcurv_tn=P_SecMean_tn(locs(1:Npeaks),:);

%Test
Fig4=figure;
Ax4=axes(Fig4);
hold(Ax4,'on'); grid(Ax4,'on');
plot(Curvature);
plot(Ax4.XLim,[Curv_Threshold,Curv_Threshold],'r','linewidth',2,'linestyle','--');
scatter(Ax4,locs,pks,100,'k','marker','*')

%Test
Fig5=figure;
Ax5=axes(Fig5);
hold(Ax5,'on'); grid(Ax5,'on');
ColorPlot(Ax5,P_SecMean_tn(:,1),P_SecMean_tn(:,2),Curvature_RGB);
scatter(Ax5,Pcurv_tn(:,1),Pcurv_tn(:,2),200,'d','marker','d')
scatter(P_tn(:,1),P_tn(:,2),0.5,'r','filled','parent',Ax5);
%% Create point vector and plot on worm in orginial orientation
[~,Pcurv_Ind]=sort(Pcurv_tn(:,1));
WormPoints_tn=[ P_SecMean_tn(1,:);Pcurv_tn(Pcurv_Ind,:);P_SecMean_tn(end,:)];
WormPoints=(coeff*WormPoints_tn')'+CG;

%Test
Fig6=figure;
Ax6=axes(Fig6);
hold(Ax6,'on');
imshow(Worm,'parent',Ax6);
plot(Ax6,WormPoints(:,1),WormPoints(:,2),'linewidth',2,'marker','d','color','r')
%% functions 
function q2=ProjectPointsOnLine(p,q1,t)
%Projects point p on line which pases through q1 with normalized direction
%t to create points q2.
%points are represented [x,y] and t the tangent normalized direction is
%[x,y] aswell.

r=p-q1;
q2=q1+r.*t;
end
function K=NumCurvature(x,y)
%Calculate curvature of an x,y curve numerically.
%returns numeric array of length(x) where K(i) is the curvature of point
%[x(i),y(i)]

n=length(x);
K=zeros(n,1);
for i=2:n-1
   x1=x(i-1);   x2=x(i);    x3=x(i+1);
   y1=y(i-1);   y2=y(i);    y3=y(i+1);
   K(i)=P3Curvature(x1,x2,x3,y1,y2,y3);
end
%assume constant curvature at ends
if n>1
    K(1)=K(2);
    K(end)=K(end-1);
end
end
function KP=P3Curvature(x1,x2,x3,y1,y2,y3)
%Calculates curvature of a point given two other adjcent points (numeric)
%https://www.mathworks.com/matlabcentral/answers/284245-matlab-code-for-computing-curvature-equation
   a = sqrt((x1-x2)^2+(y1-y2)^2); % The three sides
   b = sqrt((x2-x3)^2+(y2-y3)^2);
   c = sqrt((x3-x1)^2+(y3-y1)^2);
   A = 1/2*abs((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)); % Area of triangle
   KP = 4*A/(a*b*c); % Curvature of circumscribing circle
end
function [Tang,Norm]=NormTangLine(X,Y,k,L)
%X,Y are cordiantes. k is an index of Yk,Xk point.
%L is the length of the normal and tanget lines to be output

%output:
%Tang - [X,Y] matrix to be used with plot(Tang(:,1),Tang(:,2)
%Norm - same as tang but for normal line

dY=gradient(Y,X);
TangSlope=dY(k);

XTangDist=sqrt(L^2/(1+TangSlope^2));
XTang=[XTangDist+X(k),-XTangDist+X(k)];
YTang=(XTang-X(k))*TangSlope+Y(k);
Tang=[XTang(:),YTang(:)];

NormSlope=-1/TangSlope;
XNormDist=sqrt(L^2/(1+NormSlope^2));
XNorm=[XNormDist+X(k),-XNormDist+X(k)];
YNorm=(XNorm-X(k))*NormSlope+Y(k);
Norm=[XNorm(:),YNorm(:)];
end
function RGB=Curvature2RGB(K)
%returns [R-G-B] KX3 matrix corresponding to values of K
cmap = jet;
% make it into a index image.
cmin = min(K);
cmax = max(K);
m = length(cmap);
index = fix((K-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB=ind2rgb(index,cmap);
RGB=reshape(RGB,length(K),3);
RGB=RGB/max(RGB(:));

if std(K)<1000*eps %range of computational error - for circle and lines
    RGB=mean(RGB).*ones(size(RGB)); %turn the whole color matrix to a single color
end
end
function h=ColorPlot(axes,X,Y,RGB)
%X,Y are vectors of length n. RGB is a matrix nx3 representing color
%plots lines colored lines and returns handle array
n=length(X);
if n==1 %single point to be drawn
    scatter(axes,X,Y,100,'r')
else
    for i=1:n-1 %multiple points to be drawn as lines
        C=mean([RGB(i,:);RGB(i,:)],1);
        h(i)=plot(axes,[X(i),X(i+1)],[Y(i),Y(i+1)],'color',C,'linewidth',2);
    end
end
end