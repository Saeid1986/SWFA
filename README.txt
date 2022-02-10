%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      SPATIALLY-WEIGHTED      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                      FACTOR                        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                 ANALYSIS                  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 This is the main m-file that implements the spatially-weighted factor
 analysis (SWFA) according to the given spatially-weighted covariance
 matrix calculated by a distance-based spatio-geological weighting
 function. A work similar in application but different in mathematics can
 be found at:
 Cheng, Qiuming, Greame Bonham-Carter, Wenlei Wang, Shengyuan Zhang,
 Wenchang Li, and Xia Qinglin. "A spatially weighted principal component
 analysis for multi-element geochemical data for mapping locations of
 felsic intrusions in the Gejiu mineral district of Yunnan, China."
 Computers & Geosciences 37, no. 5 (2011): 662-669.



 The input data/parameters include:

 data: an Excel file comprising an m * p elemental data matrix, where
 each row (i = 1,...,m) represents a geochemical sample, and each column
 (j = 1,...,p) expresses a geochemical element;

 coord: an Excel file comprising the spatial locations of the samples; it
 is a three-column (XYZ) matrix, each row expressing a geochemical sample,
 whereas X, Y, and Z denoting the longitude, latitude, and elevation of
 the geochemical samples, respectively;

 rp: an Excel file including the 3D spatial locations of the potential
 sources; it is a three-column matrix, each row expresses a reference
 point (potential source), whereas X, Y, and Z denote the longitude,
 latitude, and elevation of the potential reference points, respectively;

 wl: an Excel file containing a vector of size 'm' to represent the
 geologist-defined ore productivity scores assigned to the geochemical
 samples to account for the spatial complexities of various lithological
 units within the study area;

 wh: an Excel file containing a vector of size 'm' to represent the
 geologist-defined ore productivity scores assigned to the geochemical
 samples to account for the spatial complexities of various hydrothermal
 alterations zones within the study area;

 a: an user-defined specific value to tune the relative importance of
 lithology scores;

 b: an user-defined specific value to tune the relative importance of
 hydrothermal alteration scores;

 eps: an user-specified small constant in meters subject to eps > 1 in
 order to avoid the ambiguity of the denominator;

 marg: an user-specified value in meters to decorate the margins of plots.



 The program makes it possible to apply two well-established logratio
 transformations, i.e., centered logratio (clr) and isometric logratio
 (ilr), in order to address the closure effect of the geochemical data. It
 is followed by a standardization with z-transform, which results in a
 dataset with zero means and unit variances to avoid incomparable
 magnitudes.



 The program provides the following outputs/results:

 figure 1: a plot indicating the templete of shortest distance factor;

 figure 2: a plot indicating the lithological weighting templete;

 figure 3: a plot indicating the alteration weighting templete;

 figure 4: a plot indicating the templete of overall spatial weighting
 function;

 W: an Excel file including a vector of size 'm' which represents the
 final spatial weights;

 Spatially-Weighted Covariance Matrix: an Excel file containing the
 calculated spatially-weighted covariance matrix of size p * p (for
 clr-transformed data) or p-1 * p-1 (for ilr-transformed data);

 Standard Covariance Matrix: an Excel file containing the calculated
 standard covariance matrix of size p * p (for clr-transformed data) or
 p-1 * p-1 (for ilr-transformed data);

 The program then performs the spatially-weighted or standard factor
 analysis based on the covariance matrix of interest (spatially-weighted
 or standard). This m-file deals with the principal component solution of
 the factor model through the covariance matrix (without the matrix of
 data), the latent root criterion, and uses the varimax factor rotation.
 The main reference is:
 Rencher, A. C. "Methods of Multivariate Analysis." 2nd. edition,
 New-Jersey: John Wiley & Sons. Chapter 13, (2002): 408-450.
 The program asks user-friendly many questions and returns the following
 outputs based on the user's answers:

 Complete factor analysis results such as:
 - Table of general structure of extracted components;
 - Table of unrotated principal components of the factor analysis;
 - Table of proportion of total (standardized) sample variance;
 - Table of cumulative proportion of total (standardized) sample variance;
 
 as well as optional results such as:
 - Table of varimax rotated principal components of the factor analysis;
 - Residual matrix.



 To cite this program, this would be an appropriate format:
 Esmaeiloghli, Saeid, Seyed Hassan Tabatabaei, Emmanuel John M. Carranza,
 Shahram Hosseini, and Yannick Deville. "Spatially-Weighted Factor
 Analysis for Extraction of Source-Oriented Mineralization Feature in 3D
 Coordinates of Surface Geochemical Signal." Natural Resources Research
 30, no. 6 (2021): 3925-3953.



Saeid ESMAEILOGHLI
February 26, 2019
 ------------------------
 Institut de Recherche en Astrophysique et Planétologie (IRAP),
 Centre National de la Recherche Scientifique (CNRS),
 Centre National d'Etudes Spatiales (CNES),
 Observatoire Midi-Pyrénées (OMP),
 Université Paul Sabatier (UPS),
 14 Avenue Edouard Belin,
 31400 Toulouse,
 France.
 ------------------------
 Department of Mining Engineering,
 Isfahan University of Technology (IUT),
 University Boulevard, Esteghlal Square,
 84156-83111 Isfahan,
 Iran.
 ------------------------
 Phone: (+33) 6 25 47 67 03
 E-mail: saeid@mail.fr
 ORCID: 0000-0002-7786-657X
 Website: http://esmaeiloghli.mining.iut.ac.ir
 GoogleScholar: http://scholar.google.com/citations?user=EZKZcwQAAAAJ&hl=en
 ResearchGate: https://www.researchgate.net/profile/Saeid_Esmaeiloghli
 ACADEMIA: https://cnrs.academia.edu/SaeidEsmaeiloghli