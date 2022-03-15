#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP
#define EIGEN_NO_DEBUG
#include "Eigen/Eigen"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//Reads a text file and converts it into a Matrix, given the file name and number of columns
MatrixXd Matrix_readIn(int col, string filename);
//Converts a matrix in [dms] to a vector in [decimal degrees]
VectorXd DMS_DD(MatrixXd m);
//Outputs a Matrix into a text file, given the name. Overloaded to include units
void MatToText(MatrixXd M, string filename);
void MatToText(MatrixXd M, string filename, vector<string> units);
//Another overload to put units on the first row
void MatToText(MatrixXd M, string filename, vector<string> units, char a);
//Converts a string array into a text file

//Class to store absolute orientation parameters
class AO_Params {
public:
    double tX, tY, tZ, Lambda, Omega, Phi, Kappa;
    MatrixXd M_m_to_o;
    VectorXd vHat, xHat;
    MatrixXd vHat_Mat;
    VectorXd RMSE;
    MatrixXd Cx, Corr;
};
class EOPs {
public:
    MatrixXd ro_PC_left, ro_PC_right;
    //Rotations matrices where the transformation is shown as M_from_to
    MatrixXd M_m_i_right, M_m_o, M_o_i_right, M_o_i_left;
    VectorXd omegaphikappa_left, omegaphikappa_right;
};

MatrixXd R1(double angle);
MatrixXd R2(double angle);
MatrixXd R3(double angle);

//Get initial x0 parameters to begin the adjustment
VectorXd Get_Initial_x0(MatrixXd rm, MatrixXd ro);
//Design matrix computations. First one computes the differentials of one set of XO,YO,ZO, and the second function implements them
MatrixXd Design_Matrix_EachPoint(VectorXd rm, VectorXd x0, MatrixXd M);
MatrixXd Design_Matrix(MatrixXd r_obj, MatrixXd r_mod, VectorXd x0, MatrixXd M);
//Computing Misclosure Vector
VectorXd Misclosure_Calc(VectorXd x0, MatrixXd M, MatrixXd rm, MatrixXd ro);
//Obtaining a weight matrix. Point precision, sigma_i, is assumed to be 10 micrometers
MatrixXd Create_Weight(MatrixXd ro, double c, double h, double H, double sigma_i, double B_H);
//Checking Criteria using a specified threshold for scale, rotation, and translation
bool xHat_within_tolerance(VectorXd delta, double trans, double rot, double scale);
//Using the above equations to obtain absolute orientation parameters. Base to height ratio is assumed to be 0.6
void Absolute_Orientation(AO_Params &Params, string ModelSpaceFile, string ObjectSpaceFile, double c, double h, double H, double sigma_i, double B_H=0.6);
//Converting Tie and Check points from model to Object space
void Transform_Points(AO_Params AOP, string Model_File, MatrixXd& Obj);
void Transform_Points(AO_Params AOP, MatrixXd Pts_model, MatrixXd& Obj);
void Get_M_o_i(MatrixXd M_m_o, EOPs& EOP);
void Trans_EOPs(AO_Params AOP, string ROP_file, EOPs& EOP);
void Transform_EOPs(AO_Params AO, EOPs& EOP);

#endif