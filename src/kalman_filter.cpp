#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  MatrixXd Ft = F_.transpose();

  x_ = F_ * x_;
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  long x_dim = x_.size();
  MatrixXd Id = MatrixXd::Identity(x_dim,x_dim);

  x_ = x_ + K * y;
  P_ = (Id - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  double rho = sqrt(px*px + py*py);

  if(fabs(px) < 0.00001){
    px = 0.00001;
  }
  if(fabs(rho) < 0.00001){
    rho = 0.00001;
  }

  double phi = atan2(py,px);
  if(py==0 && px==0){
    phi = 0;
  }
  double rho_pt = (px*vx + py*vy)/rho;

  VectorXd zp = VectorXd(3);
  zp << rho, phi, rho_pt;

  VectorXd y = z - zp;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  long x_dim = x_.size();
  MatrixXd Id = MatrixXd::Identity(x_dim,x_dim);

  x_ = x_ + K * y;
  P_ = (Id - K * H_) * P_;

}
