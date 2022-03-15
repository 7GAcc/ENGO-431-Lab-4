#include "Header.h"

int main() {
	double Test_c = 152.15;
	double Test_h = 272;
	double Test_H = 1462;
	double c = 153.358;
	double h = 1089.347;
	double H = 751.544;
	double test_sigma_i = 0.00001;
	//sigma i Obtained from lab 2:
	double sigma_i = 0.0000099;
	string Test_GCP_Model_Coords = "in_test_Model_Coords.txt";
	string Test_GCP_Object_Coords = "in_test_knownObj_Coords.txt";
	AO_Params Test_AO;
	EOPs test_EOP;
	
	

	cout << "Test Absolute Orientation: " << endl;
	Absolute_Orientation(Test_AO, Test_GCP_Model_Coords, Test_GCP_Object_Coords, Test_c, Test_h, Test_H, test_sigma_i);
	MatrixXd Test_Obj_Space;
	MatrixXd Test_PC_Obj_Space;
	Transform_Points(Test_AO, Test_GCP_Model_Coords, Test_Obj_Space);
	Transform_Points(Test_AO, "in_test_PC_Image.txt", Test_PC_Obj_Space);
	Transform_EOPs(Test_AO, test_EOP);

	cout << "Object Space Coordinates test: " << endl << Test_Obj_Space << endl;
	cout << "M (model to object): " << endl << Test_AO.M_m_to_o << endl;

	cout << "test M (object to Left image): " << endl << test_EOP.M_o_i_left << endl;
	cout << " And it's associated omega, phi, kappa values" << endl << test_EOP.omegaphikappa_left << endl;

	cout << "test M (object to Right image): " << endl << test_EOP.M_o_i_right << endl;
	cout << " And it's associated omega, phi, kappa values" << endl << test_EOP.omegaphikappa_right << endl;

	string GCP_Model_file = "in_CtrlPt_Model_Obs.txt";
	string GCP_Object_file = "in_CtrlPt_Obj_Coords.txt";
	string ROP_file = "in_RelOrientParams.txt";
	string Tie_Model_file = "in_Tie_model.txt";
	string Check_Model_file = "in_Check_model.txt";

	cout << "----------Lab 4 Absolute Orientation----------" << endl;
	AO_Params Lab4_AO;
	//Getting parameters
	Absolute_Orientation(Lab4_AO, GCP_Model_file, GCP_Object_file, c, h, H,sigma_i);

	//Applying parameters to other points
	MatrixXd Tie_Obj, Check_Obj, Control_Obj;
	Transform_Points(Lab4_AO, Tie_Model_file, Tie_Obj);
	Transform_Points(Lab4_AO, Check_Model_file, Check_Obj);
	Transform_Points(Lab4_AO, GCP_Model_file, Control_Obj);
	Control_Obj.conservativeResize(Control_Obj.rows() - 1, 3);

	cout << "Tie Points: " << endl << Tie_Obj << endl;
	cout << "Check Points: " << endl << Check_Obj << endl;
	EOPs Lab4_EOP;
	Trans_EOPs(Lab4_AO, ROP_file, Lab4_EOP);
	cout << "ro_PC_left: " << endl << Lab4_EOP.ro_PC_left << endl;
	cout << "ro_PC_right: " << endl << Lab4_EOP.ro_PC_right << endl;

	cout << "M (object to Left image): " << endl << Lab4_EOP.M_o_i_left << endl;
	cout << " And it's associated omega, phi, kappa values" << endl << Lab4_EOP.omegaphikappa_left << endl;

	cout << "M (object to Right image): " << endl << Lab4_EOP.M_o_i_right << endl;
	cout << " And it's associated omega, phi, kappa values" << endl << Lab4_EOP.omegaphikappa_right << endl;

	MatrixXd Check_diff=Matrix_readIn(3, "in_Check_ObjSpace_GIVEN.txt")-Check_Obj;

	//Printing everything to a text file
	MatToText(Lab4_AO.xHat, "out_AO_Parameters.txt");
	MatToText(Lab4_AO.Cx, "out_AO_Cx.txt");
	MatToText(Lab4_AO.Corr, "out_AO_Correlation.txt");
	MatToText(Lab4_AO.vHat_Mat, "out_AO_Residuals.txt");
	MatToText(Lab4_AO.RMSE, "out_AO_RMSE.txt");
	MatToText(Tie_Obj, "out_Tie_ObjSpace.txt");
	MatToText(Control_Obj, "out_Control_ObjSpace.txt");
	MatToText(Check_Obj, "out_Check_ObjSpace.txt");
	MatToText(Lab4_EOP.ro_PC_left, "out_PC_Left_ObjSpace.txt");
	MatToText(Lab4_EOP.ro_PC_right, "out_PC_right_ObjSpace.txt");
	MatToText(Check_diff, "out_Check_Pt_Difference.txt");
	MatToText(Lab4_EOP.M_o_i_left, "out_M_o_i_left.txt");
	MatToText(Lab4_EOP.omegaphikappa_left, "out_Moi_Angles_left.txt");
	MatToText(Lab4_EOP.M_o_i_right, "out_M_o_i_right.txt");
	MatToText(Lab4_EOP.omegaphikappa_right, "out_Moi_Angles_right.txt");


	MatToText(Test_AO.xHat, "out_test_AO_Parameters.txt");
	MatToText(Test_AO.Corr, "out_test_AO_Correlation.txt");
	MatToText(Test_AO.vHat_Mat, "out__test_AO_Residuals.txt");
	MatToText(Test_AO.RMSE, "out_test_AO_RMSE.txt");
	MatToText(test_EOP.M_o_i_left, "out_test_M_o_i_left.txt");
	MatToText(test_EOP.omegaphikappa_left, "out_test_Moi_Angles_left.txt");
	MatToText(test_EOP.M_o_i_right, "out_test_M_o_i_right.txt");
	MatToText(test_EOP.omegaphikappa_right, "out_test_Moi_Angles_right.txt");
	MatToText(Test_Obj_Space, "out_test_ObjSpace.txt");
	MatToText(Test_PC_Obj_Space, "out_test_PC.txt");

	system("pause");
	return 0;
}