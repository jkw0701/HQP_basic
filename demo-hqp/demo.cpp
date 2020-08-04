#include <iostream>
#include "vrep_bridge/vrep_bridge.h"

//for controller 
#include "controller/Inverse-dynamics.h"

// for tasks
#include "tasks/task-operational.h"
#include "tasks/task-joint-posture.h"
#include "tasks/task-joint-bounds.h"
#include "tasks/task-mobile.h"

// for trajectories 
#include "trajectories/trajectory-operationalspace.h"
#include "trajectories/trajectory-jointspace.h"

// for solver
#include "solvers/solver-HQP-factory.hxx"
#include "solvers/solver-utils.h"
//#include "solvers/solver-HQP-eiquadprog.h"
#include "solvers/solver-HQP-qpoases.h"

#include "utils/container.h"

#include <string>
#include <vector>
//#include <conio.h> // for keyboard hit

// for vrep keyboard event
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

// for reading files
#include <fstream>
#include <iostream>

bool _kbhit()
{
    termios term;
    tcgetattr(0, &term);

    termios term2 = term;
    term2.c_lflag &= ~ICANON;
    tcsetattr(0, TCSANOW, &term2);

    int byteswaiting;
    ioctl(0, FIONREAD, &byteswaiting);

    tcsetattr(0, TCSANOW, &term);

    return byteswaiting > 0;
}

HQP::robot::RobotModel * robot_;
HQP::InverseDynamics * invdyn_, *invdyn_2, *invdyn_3, *invdyn_total_,*invdyn2_,*invdyn3_,*invdyn_total2_ ,*invdyn4_ ;
HQP::tasks::TaskJointPosture * jointTask, *jointTask2;
HQP::tasks::TaskOperationalSpace * moveTask, * move2Task;
HQP::tasks::TaskJointLimit * jointLimitTask;
HQP::tasks::TaskMobile * mobileTask;

HQP::trajectories::TrajectoryJointCubic * trajPosture,*trajPosture2;
HQP::trajectories::TrajectoryJointConstant * trajPostureConstant;
HQP::trajectories::TrajectoryOperationCubic * trajEECubic, *trajmobile;
HQP::trajectories::TrajectoryOperationConstant * trajEEConstant;

VectorXd q(dof);
VectorXd qdot(dof);
VectorXd qdes(dof);
VectorXd q_lb(dof+2); // mobile 2 + robot 7
VectorXd q_ub(dof+2); // mobile 2 + robot 7
VectorXd wheel_prev(2); 
VectorXd wheel_curr(2); 
VectorXd tau(dof);
VectorXd dv;
double vrep_time = 0.0;
double Hz = 1000.0;
int na;
int nv;
int nq;
bool add_ = true;
using namespace HQP;
using namespace std;
double beta1;
double beta2;
solver::HQPData HQP_data1;
solver::HQPData HQP_data2;
solver::HQPData HQP_data3;

solver::HQPOutput sol_1;
solver::HQPOutput sol_2;
solver::HQPOutput sol_total_;


#define time_check
#ifdef time_check
#include <time.h>
#include <stdint.h>
#include <stdlib.h>

#define BILLION 1000000000L

int localpid(void) {
 static int a[9] = { 0 };
 return a[0];
}
#endif

int main()
{
  	
	wheel_prev.setZero();
	wheel_curr.setZero();

   robot_ = new HQP::robot::RobotModel(1); // 0: Manipulator, 1: Mobile Manipulaotr, 2: humanoid
   na = robot_->na();	
   nv = robot_->nv();

   invdyn_ = new HQP::InverseDynamics(*robot_);
   q.setZero();
   qdot.setZero();	

   solver::SolverHQPBase * solver = solver::SolverHQPFactory::createNewSolver(solver::SOLVER_HQP_QPOASES, "solver-eiquadprog1");

	// Level 0 : Joint Velocity Limit for Mobile + Manipulator 
	q_lb = -166.0 / 180.0 * M_PI * VectorXd(dof+2).setOnes();
	q_ub = -1.0*q_lb;

	q_lb.head(2) = -500.0 * VectorXd(2).setOnes();
	q_ub.head(2) = 500.0 * VectorXd(2).setOnes();

	q_lb(3) = -101.0/180.0*M_PI;
	q_ub(3) = 101.0/180.0*M_PI;

	q_lb(5) = -176.0/180.0*M_PI;
	q_ub(5) = 4.0/180.0*M_PI;

	q_lb(7) = -1.0/180.0*M_PI;
	q_ub(7) = 215/180.0*M_PI;

	double kp_jointlimit = 100.0, w_jointlimit = 1.00;
	jointLimitTask = new tasks::TaskJointLimit("joint_limit_task", *robot_);
	jointLimitTask->Kp(kp_jointlimit*VectorXd::Ones(robot_->nv()));
	jointLimitTask->Kd(2.0*jointLimitTask->Kp().cwiseSqrt());
	jointLimitTask->setJointLimit(q_lb, q_ub);

	////////////////// Joint Posture Task ////////////////////
	jointTask = new tasks::TaskJointPosture("joint_control_task", *robot_);
	double kp_posture = 400.0, w_posture = 1.00;
	jointTask->Kp(kp_posture*VectorXd::Ones(robot_->nv()-2));
	jointTask->Kd(2.0*jointTask->Kp().cwiseSqrt());

	////////////////// Operational Task ////////////////////
	moveTask = new tasks::TaskOperationalSpace("end_effector_task", *robot_, 7);
	double kp_move = 1000.0, w_move = 1.0;
	VectorXd a = VectorXd::Ones(6);
	//a.tail(3) *= 10.0;
	moveTask->Kp(kp_move*a);
	moveTask->Kd(2.0*moveTask->Kp().cwiseSqrt());
	moveTask->setSingular(false);

	//////////////////// Mobile Task ///////////////////////
	mobileTask = new tasks::TaskMobile("mobile_task", *robot_);
	double kp_mobile = 100.0, w_mobile = 1.0;
	mobileTask->Kp(kp_mobile*VectorXd::Ones(6));
	mobileTask->Kd(2.0*mobileTask->Kp().cwiseSqrt());
	mobileTask->setOnlyOriCTRL(false);


	trajPosture = new trajectories::TrajectoryJointCubic("joint_traj");
	trajEECubic = new trajectories::TrajectoryOperationCubic("operational_traj");

	trajectories::TrajectorySample sampleJoint(robot_->nv()-2);
	trajectories::TrajectorySample sampleJoint2(robot_->nv()-2);
	trajectories::TrajectorySample s(12, 6);
	trajectories::TrajectorySample s_mobile(12, 6);
	trajectories::TrajectorySample sampleEE(12, 6);

	// for v-rep
	VRepBridge vb;
	vb.isSimulationRun = false;
	vb.exitFlag = false;

	int ctrl_mode = 0;
	bool mode_change = true;
	bool HQP_flag = false;
	bool flag = false;

	double start_time = 0.0;
	VectorXd q_init(dof);


#ifdef time_check
	uint64_t diff;
 	struct timespec start, end;
 	int i;
#endif

	while (vb.simConnectionCheck() && !vb.exitFlag)
	{
		if (_kbhit()) {
			int key;
			key = getchar();
			switch (key)
			{
			case '\t':
				if (vb.isSimulationRun) {
					cout << "Simulation Pause" << endl;
					vb.isSimulationRun = false;
				}
				else {
					cout << "Simulation Run" << endl;
				//	vb._cntt = 0;
					vb.isSimulationRun = true;
				}
				break;
			case 'g': // for gravity compensation
				ctrl_mode = 0; 
				mode_change = true;
				HQP_flag = true;
				break;
			case 'h': // for home position joint ctrl
				ctrl_mode = 1;
				mode_change = true;
				HQP_flag = true;
				break;
			case 'q':
				vb.isSimulationRun = false;
				vb.exitFlag = true;
				simxEndDialog(vb.getClientID(), vb.dialog_handle, simx_opmode_oneshot);
				break;

			default:
				break;
			}
		}
		if (vb.isSimulationRun)
		{
			vb.read();
			
			if (vb._cntt > 0)
			{
				vrep_time = vb._cntt / Hz;
				VectorXd q_current(robot_->na() + 5), qdot_current(robot_->na() + 5);			

				q_current.setZero(); 
				q_current.head<2>() = vb.H_transform_.translation().head(2);
				q_current(2) = vb.euler_(2);
				q_current(3) = 0.0;
				q_current(4) = 0.0;			
				q_current.tail<dof>() = vb.current_q_;
				
				qdot_current.setZero();
				qdot_current.head<2>() = vb.H_vel_.linear().head(2);
				qdot_current(2) = vb.H_vel_.angular()(2);
				qdot_current(3) = vb.current_base_vel_(0);
				qdot_current(4) = vb.current_base_vel_(1);
				qdot_current.tail(dof) = vb.current_qdot_;
				
				robot_->getUpdateKinematics(q_current, qdot_current);
				robot_->getUpdateVelocity(vb.EE_vel_.linear(), vb.EE_vel_.angular());
				robot_->getMobilePos(vb.mobile_center);

				if (ctrl_mode == 0 ){
					if (mode_change){
						cout << "Gravity compensation" << endl;
						mode_change = false;
					}

					const VectorXd & tau = robot_->getNLEtorque();	
					//cout <<  robot_->getNLEtorque().transpose() << endl;
					vb.desired_torque_ = tau.tail(dof);
				}
				else if (ctrl_mode == 1){
					if (mode_change)
					{
						start_time = vrep_time;

						Transform3d T_endeffector;
						T_endeffector = robot_->getTransformation(7);
						trajEECubic->setInitSample(T_endeffector);

						T_endeffector.translation()(0) += 1.0;
						T_endeffector.translation()(1) += 1.0;
						
						trajEECubic->setGoalSample(T_endeffector);
						trajEECubic->setDuration(10.0);
						trajEECubic->setStartTime(start_time);
						trajEECubic->setReference(T_endeffector);			

						q_init = q_current.tail(dof);
						start_time = vrep_time;

						qdes.setZero();
						trajPosture = new trajectories::TrajectoryJointCubic("joint_traj");
						trajPosture->setInitSample(q_init);
						trajPosture->setGoalSample(q_init);
						trajPosture->setDuration(1.0);
						trajPosture->setStartTime(start_time);
						trajPosture->setReference(q_init);

						invdyn_ = new HQP::InverseDynamics(*robot_);
						invdyn_->addJointLimitTask(*jointLimitTask, 1.0, 0.0, 0.0); 
						invdyn_->addOperationalTask(*moveTask, 1.0, 1.0, 0.0);
						invdyn_->addJointPostureTask(*jointTask, 1.0, 2, 0.0); //weight, level, duration

						mode_change = false;
					}
			
					trajEECubic->setCurrentTime(vrep_time);
					s = trajEECubic->computeNext();
					moveTask->setReference(s);

					trajPosture->setCurrentTime(vrep_time);
					sampleJoint = trajPosture->computeNext();
					jointTask->setReference(sampleJoint);

					const solver::HQPData &HQPData = invdyn_->computeProblemData(vrep_time, q_current, qdot_current);
					const solver::HQPOutput &sol = solver->solve(HQPData);

					if (vb._cntt % 50 == 0)
					{
						cout << solver::HQPDataToString(HQPData, true) << endl;
					}

					const VectorXd &tau = invdyn_->getActuatorForces(sol);
					VectorXd dv = invdyn_->getAccelerations(sol);
					wheel_curr = dv.head(2);

					vb.desired_torque_ = tau;
					vb.desired_base_vel_(0) = wheel_curr(0);
					vb.desired_base_vel_(1) = wheel_curr(1);
					vb.desired_base_vel_(2) = wheel_curr(1);
					vb.desired_base_vel_(3) = wheel_curr(0);

					wheel_curr = wheel_prev;
				}

				vb.write();

			}
			
			vb.simLoop();
			vb._cntt++;

		}
	}

	return 0;
}
