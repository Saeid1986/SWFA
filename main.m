%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     SPATIALLY-WEIGHTED      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%            FACTOR             %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%          ANALYSIS           %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Author: Saeid ESMAEILOGHLI
% February 26, 2019
% ------------------------
% Institut de Recherche en Astrophysique et Planétologie (IRAP),
% Centre National de la Recherche Scientifique (CNRS),
% Centre National d'Etudes Spatiales (CNES),
% Observatoire Midi-Pyrénées (OMP),
% Université Paul Sabatier (UPS),
% 14 Avenue Edouard Belin,
% 31400 Toulouse,
% France.
% ------------------------
% Department of Mining Engineering,
% Isfahan University of Technology (IUT),
% University Boulevard, Esteghlal Square,
% 84156-83111 Isfahan,
% Iran.
% ------------------------
% Phone: (+33) 6 25 47 67 03
% E-mail: saeid@mail.fr
% ORCID: 0000-0002-7786-657X
% Website: http://esmaeiloghli.mining.iut.ac.ir
% GoogleScholar: http://scholar.google.com/citations?user=EZKZcwQAAAAJ&hl=en
% ResearchGate: https://www.researchgate.net/profile/Saeid_Esmaeiloghli
% ACADEMIA: https://cnrs.academia.edu/SaeidEsmaeiloghli


% -------------------------------------------------------------------------

% Author's Words:
%
% This is the main m-file that implements the spatially-weighted factor
% analysis (SWFA) according to the given spatially-weighted covariance
% matrix calculated by a distance-based spatio-geological weighting
% function. A work similar in application but different in mathematics can
% be found at:
% Cheng, Qiuming, Greame Bonham-Carter, Wenlei Wang, Shengyuan Zhang,
% Wenchang Li, and Xia Qinglin. "A spatially weighted principal component
% analysis for multi-element geochemical data for mapping locations of
% felsic intrusions in the Gejiu mineral district of Yunnan, China."
% Computers & Geosciences 37, no. 5 (2011): 662-669.
%
%
%
% The input data/parameters include:
%
% data: an Excel file comprising an m * p elemental data matrix, where
% each row (i = 1,...,m) represents a geochemical sample, and each column
% (j = 1,...,p) expresses a geochemical element;
%
% coord: an Excel file comprising the spatial locations of the samples; it
% is a three-column (XYZ) matrix, each row expressing a geochemical sample,
% whereas X, Y, and Z denoting the longitude, latitude, and elevation of
% the geochemical samples, respectively;
%
% rp: an Excel file including the 3D spatial locations of the potential
% sources; it is a three-column matrix, each row expresses a reference
% point (potential source), whereas X, Y, and Z denote the longitude,
% latitude, and elevation of the potential reference points, respectively;
%
% wl: an Excel file containing a vector of size 'm' to represent the
% geologist-defined ore productivity scores assigned to the geochemical
% samples to account for the spatial complexities of various lithological
% units within the study area;
%
% wh: an Excel file containing a vector of size 'm' to represent the
% geologist-defined ore productivity scores assigned to the geochemical
% samples to account for the spatial complexities of various hydrothermal
% alterations zones within the study area;
%
% a: an user-defined specific value to tune the relative importance of
% lithology scores;
%
% b: an user-defined specific value to tune the relative importance of
% hydrothermal alteration scores;
%
% eps: an user-specified small constant in meters subject to eps > 1 in
% order to avoid the ambiguity of the denominator;
%
% marg: an user-specified value in meters to decorate the margins of plots.
%
%
%
% The program makes it possible to apply two well-established logratio
% transformations, i.e., centered logratio (clr) and isometric logratio
% (ilr), in order to address the closure effect of the geochemical data. It
% is followed by a standardization with z-transform, which results in a
% dataset with zero means and unit variances to avoid incomparable
% magnitudes.
%
%
%
% The program provides the following outputs/results:
%
% figure 1: a plot indicating the templete of shortest distance factor;
%
% figure 2: a plot indicating the lithological weighting templete;
%
% figure 3: a plot indicating the alteration weighting templete;
%
% figure 4: a plot indicating the templete of overall spatial weighting
% function;
%
% W: an Excel file including a vector of size 'm' which represents the
% final spatial weights;
%
% Spatially-Weighted Covariance Matrix: an Excel file containing the
% calculated spatially-weighted covariance matrix of size p * p (for
% clr-transformed data) or p-1 * p-1 (for ilr-transformed data);
%
% Standard Covariance Matrix: an Excel file containing the calculated
% standard covariance matrix of size p * p (for clr-transformed data) or
% p-1 * p-1 (for ilr-transformed data);
%
% The program then performs the spatially-weighted or standard factor
% analysis based on the covariance matrix of interest (spatially-weighted
% or standard). This m-file deals with the principal component solution of
% the factor model through the covariance matrix (without the matrix of
% data), the latent root criterion, and uses the varimax factor rotation.
% The main reference is:
% Rencher, A. C. "Methods of Multivariate Analysis." 2nd. edition,
% New-Jersey: John Wiley & Sons. Chapter 13, (2002): 408-450.
% The program asks user-friendly many questions and returns the following
% outputs based on the user's answers:
%
% Complete factor analysis results such as:
% - Table of general structure of extracted components;
% - Table of unrotated principal components of the factor analysis;
% - Table of proportion of total (standardized) sample variance;
% - Table of cumulative proportion of total (standardized) sample variance;
% 
% as well as optional results such as:
% - Table of varimax rotated principal components of the factor analysis;
% - Residual matrix.
%
%
%
% To cite this program, this would be an appropriate format:
% Esmaeiloghli, Saeid, Seyed Hassan Tabatabaei, Emmanuel John M. Carranza,
% Shahram Hosseini, and Yannick Deville. "Spatially-Weighted Factor
% Analysis for Extraction of Source-Oriented Mineralization Feature in 3D
% Coordinates of Surface Geochemical Signal." Natural Resources Research
% 30, no. 6 (2021): 3925-3953.


%% Clear the Screen

clc
clear
close all

%% Load Input Data

data = xlsread('geochemical data');

coord = xlsread('sample coordinates');

rp = xlsread('reference points');

wl = xlsread('lithological weights');

wh = xlsread('hydrothermal alteration weights');

a = 1;

b = 1;

eps = 2;

marg = 100;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%        NOW        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   YOU CAN RUN THE   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%      PROGRAM      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% A. Centered Log-Ratio (CLR) or Isometric Log-Ratio (ILR) Transformation

y = clr(data);
% y = ilr(data);

%% B. Standardizing by z-Transform (Zero-Mean & Unit-Variance)

n = length(y(:,1));
p = length(y(1,:));

mu = mean(y);
sd = std(y);

for i = 1:n
    
    for j = 1:p
        
        z(i,j) = (y(i,j) - mu(j))/sd(j);
        
    end
   
end

%% C. Calculating the Spatial Weighting Factor (W)

X = coord(:,1);
Y = coord(:,2);
Z = coord(:,3);

n_dim = length(X);
n_r = length(rp);

for i = 1:n_dim
    
    for j = 1:n_r
        
        alldist(i,j) = sqrt(((rp(j,1) - X(i))^2) + ((rp(j,2) - Y(i))^2)...
            + ((rp(j,3) - Z(i))^2)); 
        
    end
    
end

for i = 1:n_dim
    
    mindist(i) = min(alldist(i,:));
    
end

for i = 1:n_dim
    
    ws(i) = 1./(log(mindist(i)) + eps);
    
end

ws = ws';
ns_w = ws.*(wl.^a).*(wh.^b);

for i = 1:n_dim
        
    s_w(i) = (ns_w(i) - min(ns_w))./(max(ns_w) - min(ns_w));
        
    w = s_w';
    
end

xlswrite('W.xls',w);

%% D. Plotting the Spatially Weighting Factor (W)

xmin = min(X) - marg;
xmax = max(X) + marg;
ymin = min(Y) - marg;
ymax = max(Y) + marg;

figure(1)                           % Distance weighting factor plot.

scatter(X,Y,50,ws,'filled');
colorbar
title('Shortest Distance','fontname','Cambria','FontSize',14)
xlabel('Easting','fontname','Cambria','FontSize',12)
ylabel('Northing','fontname','Cambria','FontSize',12)
xlim([xmin xmax])
ylim([ymin ymax])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(2)                           % Lithological weighting factor plot.

scatter(X,Y,50,wl,'filled');
colorbar
title('Lithology Weight','fontname','Cambria','FontSize',14)
xlabel('Easting','fontname','Cambria','FontSize',12)
ylabel('Northing','fontname','Cambria','FontSize',12)
xlim([xmin xmax])
ylim([ymin ymax])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(3)                           % Alteration weighting factor plot.

scatter(X,Y,50,wh,'filled');
colorbar
title('Hydrothermal Alteration Weight','fontname','Cambria','FontSize',14)
xlabel('Easting','fontname','Cambria','FontSize',12)
ylabel('Northing','fontname','Cambria','FontSize',12)
xlim([xmin xmax])
ylim([ymin ymax])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(4)                           % Spatial weighting factor (W) plot.

scatter(X,Y,50,w,'filled');
colorbar
title('Spatial Weighting Function','fontname','Cambria','FontSize',14)
xlabel('Easting','fontname','Cambria','FontSize',12)
ylabel('Northing','fontname','Cambria','FontSize',12)
xlim([xmin xmax])
ylim([ymin ymax])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

%% E. Calculating the Spatially-Weighted Covariance Matrix

% E-1. Calculating the weighted mean value for each variable

j_dim = length(z(1,:));
i_dim = length(z(:,1));

for j = 1:j_dim

    for i = 1:i_dim
       
        a(i,j) = w(i).*z(i,j);
        nu(j) = sum(a(:,j))/sum(w); % The weighted mean value of variables.
    
    end
   
end

% E-2. Calculating the weighted corvariance matrix

for j = 1:j_dim
       
    for k = 1:j_dim
            
        c(j,k) = sum(w(:,1).*(z(:,j) - nu(j)).*(z(:,k) - nu(k)));
        d(j,k) = sum(w(:,1));
            
        s(j,k) = c(j,k)./d(j,k);
            
    end     
    
end

xlswrite('Spatially-Weighted Covariance Matrix',s);

%% F. Calculating the Standard Covariance Matrix

C = cov(z);
xlswrite('Standard Covariance Matrix',C);

%% G. Spatially-Weighted (SWFA) & Standard Factor Analysis (FA)

fabycov(s)
