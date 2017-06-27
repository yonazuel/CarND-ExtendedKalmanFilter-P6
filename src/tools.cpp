#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse = VectorXd(4);
	rmse << 0,0,0,0;

	for (int i = 0; i<estimations.size(); ++i){
		VectorXd res = estimations[i] - ground_truth[i];
		VectorXd se = res.array() * res.array();
		rmse += se;
	}

	rmse /= estimations.size();
	return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj = MatrixXd(3,4);
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double rho2 = px*px + py*py;
	if(fabs(rho2) < 0.00001){
		rho2 = 0.00001;
	}

	double rho = sqrt(rho2);
	double rho_3_2 = rho * rho2;
	double det = vx*py-vy*px;

	Hj << px/rho, py/rho, 0, 0,
		-py/rho2, px/rho2, 0, 0,
		py*det/rho_3_2, -px*det/rho_3_2, px/rho, py/rho;

	return Hj;
}
