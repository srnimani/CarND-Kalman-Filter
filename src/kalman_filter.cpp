#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}


void KalmanFilter::Predict()
{
    /**
     TODO: Done
     * predict the state
     */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}


void KalmanFilter::Update(const VectorXd &z)
{
    /**
     TODO: Done.. Chapter 6 Lidar & Radar Fusion with KF in C++
     * Update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
    
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    
    // New estimated state
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
  TODO: Done, see Course notes from chapter 16 to 20
    * update the state by using Extended Kalman Filter equations
  */
    
    double px = x_[0];
    double py = x_[1];
    double vx = x_[2];
    double vy = x_[3];
    
    double rho = sqrt(px * px + py * py);
    double phi ;
    double rho_dot ;
    
    // Check for divide by 0..
    
    if (px != 0)
    {
        phi= atan2(py, px);
        rho_dot = (px * vx + py * vy) / rho;
    }
    
    else
    {
        phi= 0;
        rho_dot = 0;
    }
    
    VectorXd z_pred(3);
    z_pred << rho, phi, rho_dot;
    
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
    
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    
    // New state estimate
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
    
}


