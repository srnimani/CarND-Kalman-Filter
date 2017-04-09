
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Acceleration noise component

int noise_ax = 9;
int noise_ay = 9;


/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
    is_initialized_ = false;
    previous_timestamp_ = 0.0;
    
    // Define and initialize the matrices
    
    // Measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
    0, 0.0225;
    
    // Measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
    
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
    
    // State covariance matrix P
    MatrixXd P_ = MatrixXd(4, 4);
    P_ <<   1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;
    
    MatrixXd F_ = MatrixXd(4, 4);
    MatrixXd Q_ = MatrixXd(4, 4);
    
    // Vector state, unknown at the beginning
    VectorXd x_ = VectorXd(4);
    
    /**
     TODO: Done
     * Finish initializing the FusionEKF
     * Set the process and measurement noises */

    ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
    
   }

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    
    if (!is_initialized_) {
        
        /**
         TODO - Done
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */

        // first measurement, start with unknown state and all initialized to 0
        float px = 0;
        float py = 0;
        float vx = 0;
        float vy = 0;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho = measurement_pack.raw_measurements_[0]; // range value
            float phi = measurement_pack.raw_measurements_[1]; // angle between rho and x
            float rho_dot = measurement_pack.raw_measurements_[2]; // range rate
            
            px = rho * cos(phi);
            py = rho * sin(phi);
            vx = rho_dot * cos(phi);
            vy = rho_dot * sin(phi);
            
            if(fabs(px) < 0.0001)
            {
                px = 1;
                ekf_.P_(0,0) = 1000;
            }
            
            if(fabs(py) < 0.0001)
            {
                py = 1;
                ekf_.P_(1,1) = 1000;
            }
        }
        
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
            /**
             Initialize state for Laser measurements
             */
            px = measurement_pack.raw_measurements_[0];
            py = measurement_pack.raw_measurements_[1];
        }
        
        ekf_.x_ << px, py, vx, vy;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        previous_timestamp_ = measurement_pack.timestamp_;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    /**
     TODO: Done
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    
    ekf_.F_ <<  1, 0, dt, 0,
    0, 1, 0, dt,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
    // Calculate frequently used factors
    double dt_2 = pow(dt, 2);
    double dt_3 = pow(dt, 3) / 2;
    double dt_4 = pow(dt, 4) / 4;
    
    ekf_.Q_ <<  dt_4 * noise_ax, 0, dt_3 * noise_ax, 0,
    0, dt_4 * noise_ay, 0, dt_3 * noise_ay,
    dt_3 * noise_ax, 0, dt_2 * noise_ax, 0,
    0, dt_3 * noise_ay, 0, dt_2 * noise_ay;
    
    ekf_.Predict();
    
    // Update timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
   
    /**
     TODO: Done
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Radar updates and look for divide by zero error!
        try
        {
            ekf_.R_ = R_radar_;
            ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
            ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        } catch(...){
            return; // Jacobian error, ignore update
        }
        
    }
    else
    {
        // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
