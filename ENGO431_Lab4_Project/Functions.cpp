#include "Header.h"

double PI = atan(1)*4;

MatrixXd Matrix_readIn(int col, string filename) {
    ifstream infile;
    MatrixXd result(1, col);
    infile.open(filename, ifstream::in);
    int rows = 1;
    while (!infile.eof()) {
        double a;
        result.conservativeResize(rows, col);
        for (int i = 0; i < col; i++) {
            infile >> a;
            result(rows - 1, i) = a;
        }
        rows++;
    }
    infile.close();
    return result;
}

VectorXd DMS_DD(MatrixXd m) {
    VectorXd result(m.rows());
    for (int i = 0; i < m.rows(); i++) {
        result(i) = m(i, 0) + m(i, 1) / 60 + m(i, 2) / 3600;
    }
    return result;
}

//Outputs Matrix to text file
void MatToText(MatrixXd M, string filename) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << val << "\t"; //Stores i-th row and j-th column in text file
        }
        outfile << endl;
    }
    outfile.close();
}
void MatToText(MatrixXd M, string filename, vector<string> units) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << val << "\t"; //Stores i-th row and j-th column in text file
        }
        outfile << units[i] << endl;
    }
    outfile.close();
}
void MatToText(MatrixXd M, string filename, vector<string> units, char a) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.cols(); i++) {
        outfile << units[i] << "\t";
    }
    outfile << endl;
    for (int i = 0; i < M.rows(); i++) {
        if (i != 0)
            outfile << endl;
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << fixed << setprecision(9) << val << "\t"; //Stores i-th row and j-th column in text file
        }

    }
    outfile.close();
}
//Rotation in the X axis
MatrixXd R1(double angle) {
    MatrixXd result(3, 3);
	double rad = angle;// *PI / 180;
    result << 1, 0, 0,
        0, cos(rad), sin(rad),
        0, -sin(rad), cos(rad);
    return result;
}
//Rotation in the Y axis
MatrixXd R2(double angle) {
    MatrixXd result(3, 3);
	double rad = angle;// *PI / 180;
    result << cos(rad), 0, -sin(rad),
        0, 1, 0,
        sin(rad), 0, cos(rad);
    return result;
}
//Rotation in the Z axis
MatrixXd R3(double angle) {
    MatrixXd result(3, 3);
	double rad = angle;// *PI / 180;
    result << cos(rad), sin(rad), 0,
        -sin(rad), cos(rad), 0,
        0, 0, 1;
    return result;
}
//Getting initial x0 using given model and object coordinates of control points
VectorXd Get_Initial_x0(MatrixXd rm, MatrixXd ro) {
    VectorXd result(7);
    double Lambda, Omega, Phi, Kappa;
    VectorXd t;
    //In the x0 vector, the variables are in the order: Omega,Phi,Kappa,tX,tY,tZ,Lambda
    //Assuming straight level flight, therefore Omega and Phi are 0
    Omega = 0;
    Phi = 0;
    
    //Calculating Kappa using first and third control points
    //Kappa = atan2(ro(2, 0) - ro(0, 0), ro(2, 1) - ro(0, 1)) * 180 / PI - atan2(rm(2, 0) - rm(0, 0), rm(2, 1) - rm(0, 1)) * 180 / PI;
	Kappa = atan2(ro(2, 0) - ro(0, 0), ro(2, 1) - ro(0, 1)) - atan2(rm(2, 0) - rm(0, 0), rm(2, 1) - rm(0, 1));
    //Calculating Lamda using first and third control points
        //Not including the Z coordinate in this case because it is not necessary
    double do_ij = sqrt(pow(ro(2, 0) - ro(0, 0), 2) + pow(ro(2, 1) - ro(0, 1), 2) + pow(ro(2, 2) - ro(0, 2), 2));
    double dm_ij = sqrt(pow(rm(2, 0) - rm(0, 0), 2) + pow(rm(2, 1) - rm(0, 1), 2) + pow(rm(2, 2) - rm(0, 2), 2));
    Lambda = do_ij / dm_ij;
    //Calculating tX, tY, tZ using second control point

    t = ro.row(1) - Lambda * R3(Kappa) * rm.row(1);

    result << Omega, Phi, Kappa, t(0), t(1), t(2), Lambda;
    return result;
}
MatrixXd Design_Matrix_EachPoint(VectorXd rm, VectorXd x0, MatrixXd M) {
    MatrixXd result(3, 7);
    double x, y, z, omega, phi, kappa, lambda;
    x = rm(0); y = rm(1); z = rm(2);
    //omega = x0(0) * PI / 180; phi = x0(1) * PI / 180; kappa = x0(2) * PI / 180; lambda = x0(6);
	omega = x0(0); phi = x0(1); kappa = x0(2); lambda = x0(6);

    VectorXd dX(7);
    VectorXd dY(7);
    VectorXd dZ(7);
    dX(0) = lambda * y * (-sin(omega) * sin(kappa) + cos(omega) * sin(phi) * cos(kappa)) + lambda * z * (cos(omega) * sin(kappa) + sin(omega) * sin(phi) * cos(kappa));
    dX(1) = -lambda * x * sin(phi) * cos(kappa) + lambda * y * sin(omega) * cos(phi) * cos(kappa) - lambda * z * cos(omega) * cos(phi) * cos(kappa);
    dX(2) = -lambda * x * cos(phi) * sin(kappa) + lambda * y * (cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa)) + lambda * z * (sin(omega) * cos(kappa) + cos(omega) * sin(phi) * sin(kappa));
    dX(3) = 1;
    dX(4) = 0;
    dX(5) = 0;
    dX(6) = x * M(0, 0) + y * M(0, 1) + z * M(0, 2);
 
    dY(0) = lambda * y * (-sin(omega) * cos(kappa) - cos(omega) * sin(phi) * sin(kappa)) + lambda * z * (cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa));
    dY(1) = lambda * x * sin(phi) * sin(kappa) - lambda * y * sin(omega) * cos(phi) * sin(kappa) + lambda * z * cos(omega) * cos(phi) * sin(kappa);
    dY(2) = -lambda * x * cos(phi) * cos(kappa) + lambda * y * (-cos(omega) * sin(kappa) - sin(omega) * sin(phi) * cos(kappa)) + lambda * z * (-sin(omega) * sin(kappa) + cos(omega) * sin(phi) * cos(kappa));
    dY(3) = 0;
    dY(4) = 1;
    dY(5) = 0;
    dY(6) = x * M(1, 0) + y * M(1, 1) + z * M(1, 2);

    dZ(0) = -lambda * y * cos(omega) * cos(phi) - lambda * z * sin(omega) * cos(phi);
    dZ(1) = lambda * x * cos(phi) + lambda * y * sin(omega) * sin(phi) - lambda * z * cos(omega) * sin(phi);
    dZ(2) = 0;
    dZ(3) = 0;
    dZ(4) = 0;
    dZ(5) = 1;
    dZ(6) = x * M(2, 0) + y * M(2, 1) + z * M(2, 2);

    result.row(0) = dX;
    result.row(1) = dY;
    result.row(2) = dZ;

    return result;
}
MatrixXd Design_Matrix(MatrixXd r_obj, MatrixXd r_mod, VectorXd x0, MatrixXd M) {
    MatrixXd result(r_obj.rows() * 3, 7);
    MatrixXd dXdYdZ(3, 7);
    for (int i = 0; i < r_obj.rows(); i++) {
        dXdYdZ = Design_Matrix_EachPoint(r_mod.row(i), x0, M);
        result.row(3 * i) = dXdYdZ.row(0);
        result.row(3 * i + 1) = dXdYdZ.row(1);
        result.row(3 * i + 2) = dXdYdZ.row(2);
    }
    return result;
}
VectorXd Misclosure_Calc(VectorXd x0, MatrixXd M, MatrixXd rm, MatrixXd ro) {
    VectorXd result(rm.rows() * 3);
    for (int i = 0; i < rm.rows(); i++) {
        result(3 * i) = x0(6) * (M(0, 0) * rm(i, 0) + M(0, 1) * rm(i, 1) + M(0, 2) * rm(i, 2)) + x0(3) - ro(i, 0);
        result(3 * i + 1) = x0(6) * (M(1, 0) * rm(i, 0) + M(1, 1) * rm(i, 1) + M(1, 2) * rm(i, 2)) + x0(4) - ro(i, 1);
        result(3 * i + 2) = x0(6) * (M(2, 0) * rm(i, 0) + M(2, 1) * rm(i, 1) + M(2, 2) * rm(i, 2)) + x0(5) - ro(i, 2);
    }
    return result;
}
MatrixXd Create_Weight(MatrixXd ro, double c, double h, double H, double sigma_i, double B_H) {
    MatrixXd result(ro.rows() * 3, ro.rows() * 3);
    result.fill(0);

    //Getting scale
    double S = abs(H - h) / c;
    //x and y precision are assumed to be the same
    double sigma_xy = S * sigma_i;
    double sigma_z = sqrt(2) * S * sigma_i/ B_H;
    //Using the above sigmas to obtain the weight matrix
    cout << "Expected Precision in x and y (assumed to be equal): " << sigma_xy << " m" << endl;
    cout << "Expected Precision in z: " << sigma_z << " m" << endl;
    for (int i = 0; i < ro.rows(); i++) {
        result(3 * i, 3 * i) = 1 / sigma_xy;
        result(3 * i + 1, 3 * i + 1) = 1 / sigma_xy;
        result(3 * i + 2, 3 * i + 2) = 1 / sigma_z;
    }
    return result;
}
bool xHat_within_tolerance(VectorXd delta, double trans, double rot, double scale) {
    if (delta(0) < rot && delta(1) < rot && delta(2) < rot && delta(3) < trans && delta(4) < trans && delta(5) < trans && delta(6) < scale) {
        return true;
    }
    else {
        return false;
    }
}
void Absolute_Orientation(AO_Params& Params, string ModelSpaceFile, string ObjectSpaceFile, double c, double h, double H, double sigma_i, double B_H) {
//cout << "Works" << endl;
    //Reading in Model space (Constants) and object space (Observations) coordinates of GCPs from text file:
    MatrixXd rm_GCP = Matrix_readIn(3, ModelSpaceFile);
    MatrixXd ro_GCP = Matrix_readIn(3, ObjectSpaceFile);
    //Getting number of observations, unknowns, and degrees of freedom
        //m=n
    int n = rm_GCP.rows() * 3;
    int u = 7;
    int r = n - u;
    int count = 0;
    //Getting initial approximation of x0
    VectorXd x0 = Get_Initial_x0(rm_GCP, ro_GCP);
    //x0 << -0.824127, -0.717738, 18.891137, 6349.551, 3964.645, 1458.114, 7.585632;
    //cout << x0 << endl;
    MatrixXd A(n, u);
    MatrixXd M(3, 3);
    VectorXd w(n);
    VectorXd xHat(u);
    VectorXd delta(u);
    //Creating weight matrix
    MatrixXd P = Create_Weight(ro_GCP, c/1000, h, H, sigma_i, B_H);

    //Beginning the adjustment
    do {
        M = R3(x0(2)) * R2(x0(1)) * R1(x0(0));
        A = Design_Matrix(ro_GCP, rm_GCP, x0, M);
        w = Misclosure_Calc(x0, M, rm_GCP, ro_GCP);
        delta = -(A.transpose() * P * A).inverse() * A.transpose() * P * w;
        xHat = x0 + delta;
        x0 = xHat;

        //if (count > 200) {
        //    std::cout << "Threshold is too small" << endl << endl;
        //    break;
        //}

        count++;
        
    } while (xHat_within_tolerance(delta, 0.001, 1 / 3600, 0.000001) == false);

	//CONVERT OMEGA PHI KAPPA TO DEGREES
	xHat(0) = xHat(0) * 180 / PI;
	xHat(1) = xHat(1) * 180 / PI;
	xHat(2) = xHat(2) * 180 / PI;

	std::cout << "xHat after " << count << " iterations" << endl << xHat << endl;
    Params.Omega = xHat(0);
    Params.Phi = xHat(1);
    Params.Kappa = xHat(2);
    Params.tX = xHat(3);
    Params.tY = xHat(4);
    Params.tZ = xHat(5);
    Params.Lambda = xHat(6);
    Params.M_m_to_o = M;
    Params.vHat = A * delta + w;
    MatrixXd vHat(A.rows()/3, 3);
    for (int i = 0; i < vHat.rows(); i++) {
        vHat(i, 0) = Params.vHat(3 * i);
        vHat(i, 1) = Params.vHat(3 * i + 1);
        vHat(i, 2) = Params.vHat(3 * i + 2);
    }
    Params.vHat_Mat = vHat;
    VectorXd RMSE(3);
    RMSE.fill(0);
    for (int i = 0; i < vHat.rows(); i++) {
        RMSE(0) = RMSE(0) + pow(vHat.col(0).sum() / vHat.rows() - vHat(i, 0), 2);
        RMSE(1) = RMSE(1) + pow(vHat.col(1).sum() / vHat.rows() - vHat(i, 1), 2);
        RMSE(2) = RMSE(2) + pow(vHat.col(2).sum() / vHat.rows() - vHat(i, 2), 2);
    }
    RMSE /= vHat.rows();
    for (int i = 0; i < RMSE.rows(); i++) {
        RMSE(i) = sqrt(RMSE(i));
    }
    cout << "RMSE: " << endl << RMSE << endl;
    //Getting the Cx matrix and correlation matrix
    Params.Cx = (A.transpose() * P * A).inverse();
    MatrixXd Corr(Params.Cx.rows(), Params.Cx.rows());
    for (int i = 0; i < Corr.rows(); i++) {
        for (int j = 0; j < Corr.cols(); j++) {
            Corr(i, j) = Params.Cx(i, j) / sqrt(Params.Cx(i, i) * Params.Cx(j, j));
        }
    }
    Params.Corr = Corr;
    Params.xHat = xHat;
    Params.RMSE = RMSE;
}


//Second Part of AO:
void Transform_Points(AO_Params AOP, string Model_File, MatrixXd &Obj) {
    MatrixXd Pts_model = Matrix_readIn(3, Model_File);
    Obj.resize(Pts_model.rows(), 3);
    VectorXd t(3);
    t << AOP.tX, AOP.tY, AOP.tZ;
    VectorXd Vec, r_m;
    double Lam = AOP.Lambda;
    for (int i = 0; i < Obj.rows(); i++) {
        r_m = Pts_model.row(i);
        Vec = Lam * AOP.M_m_to_o * r_m + t;
        Obj.row(i) = Vec;
    }

}
void Transform_Points(AO_Params AOP, MatrixXd Pts_model, MatrixXd& Obj) {
    Obj.resize(Pts_model.rows(), 3);
    VectorXd t(3);
    t << AOP.tX, AOP.tY, AOP.tZ;
    VectorXd Vec, r_m;
    double Lam = AOP.Lambda;
    for (int i = 0; i < Obj.rows(); i++) {
        r_m = Pts_model.row(i);
        Vec = Lam * AOP.M_m_to_o * r_m + t;
        Obj.row(i) = Vec;
    }

}
void Get_M_o_i(MatrixXd M_m_o, EOPs& EOP) {
    EOP.M_o_i_right = EOP.M_m_i_right * M_m_o.transpose();
    EOP.M_o_i_left = M_m_o.transpose();
    EOP.omegaphikappa_left.resize(3);
    EOP.omegaphikappa_right.resize(3);
    double w, phi, kap;
    w = atan2(-EOP.M_o_i_left(2, 1), EOP.M_o_i_left(2, 2))*180/PI;
    phi = asin(EOP.M_o_i_left(2, 0));
    kap = atan2(-EOP.M_o_i_left(1, 0), EOP.M_o_i_left(0, 0))*180/PI;
    EOP.omegaphikappa_left << w, phi, kap;
    w = atan2(-EOP.M_o_i_right(2, 1), EOP.M_o_i_right(2, 2)) * 180 / PI;
    phi = asin(EOP.M_o_i_right(2, 0));
    kap = atan2(-EOP.M_o_i_right(1, 0), EOP.M_o_i_right(0, 0)) * 180 / PI;
    EOP.omegaphikappa_right << w, phi, kap;
}
void Trans_EOPs(AO_Params AOP, string ROP_file, EOPs& EOP) {
    MatrixXd ROP = Matrix_readIn(1, ROP_file);
    MatrixXd B(1, 3);
    B << ROP(0), ROP(1), ROP(2);
    MatrixXd rm_PC_left(1, 3);
    rm_PC_left << 0, 0, 0;
    Transform_Points(AOP, rm_PC_left, EOP.ro_PC_left);
    Transform_Points(AOP, B, EOP.ro_PC_right);
    EOP.M_m_i_right=R3(ROP(5))*R2(ROP(4))*R1(ROP(3));

    Get_M_o_i(AOP.M_m_to_o, EOP);

}
void Transform_EOPs(AO_Params AO, EOPs& EOP) {
    MatrixXd Mmo = AO.M_m_to_o;
    MatrixXd Mmi_l(3, 3);
    Mmi_l << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    MatrixXd Mmi_r(3, 3);
    Mmi_r << 0.9981, 0.0553, -0.0259,
        -0.0551, 0.9984, 0.0091,
        0.0263, -0.0077, 0.9996;
    EOP.M_o_i_right = Mmi_r * Mmo.transpose();
    EOP.M_o_i_left = Mmo.transpose();

    VectorXd OPK(3);
    OPK(0) = atan2(-EOP.M_o_i_left(2, 1), EOP.M_o_i_left(2, 2)) * 180 / PI;
    OPK(1) = asin(EOP.M_o_i_left(2, 0)) * 180 / PI;
    OPK(2)= atan2(-EOP.M_o_i_left(1, 0), EOP.M_o_i_left(0, 0)) * 180 / PI;

    EOP.omegaphikappa_left = OPK;
    
    OPK(0) = atan2(-EOP.M_o_i_right(2, 1), EOP.M_o_i_right(2, 2)) * 180 / PI;
    OPK(1) = asin(EOP.M_o_i_right(2, 0)) * 180 / PI;
    OPK(2) = atan2(-EOP.M_o_i_right(1, 0), EOP.M_o_i_right(0, 0)) * 180 / PI;

    EOP.omegaphikappa_right = OPK;

}




