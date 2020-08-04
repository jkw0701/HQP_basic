#include "trajectories/trajectory-operationalspace.h"
#include "utils/utils.h"
#include <iostream>
using namespace HQP::trajectories;
using namespace std;

TrajectoryOperationConstant::TrajectoryOperationConstant(const std::string & name)
    :TrajectoryBase(name)
{}

TrajectoryOperationConstant::TrajectoryOperationConstant(const std::string & name, const Transform3d & M) :TrajectoryBase(name)
{
    m_sample.resize(12, 6);
	m_sample.pos.head<3>() = M.translation();
	typedef Eigen::Matrix<double, 9, 1> Vector9;
	m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&M.rotation()(0), 9);
}
void TrajectoryOperationConstant::setReference(const Transform3d ref) {
	m_sample.resize(12, 6);
	m_sample.pos.head<3>() = ref.translation();
	typedef Eigen::Matrix<double, 9, 1> Vector9;
	m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&ref.rotation()(0), 9);
}
unsigned int TrajectoryOperationConstant::size() const
{
    return 6;
}

const TrajectorySample & TrajectoryOperationConstant::operator()(double )
{
    return m_sample;
}

const TrajectorySample & TrajectoryOperationConstant::computeNext()
{
    return m_sample;
}

void TrajectoryOperationConstant::getLastSample(TrajectorySample & sample) const
{
    sample = m_sample;
}

bool TrajectoryOperationConstant::has_trajectory_ended() const
{
    return true;
}
/////////////////////////////////////////////////
TrajectoryOperationCubic::TrajectoryOperationCubic(const std::string & name)
	:TrajectoryBase(name)
{}

TrajectoryOperationCubic::TrajectoryOperationCubic(const std::string & name, const Transform3d & init_M, const Transform3d & goal_M, const double & duration, const double & stime) : TrajectoryBase(name)
{
	setGoalSample(goal_M);
	setInitSample(init_M);
	setDuration(duration);
	setStartTime(stime);

	//m_sample.resize(12, 6);
	//m_sample.pos.head<3>() = init_M.translation();
	//typedef Eigen::Matrix<double, 9, 1> Vector9;
	//m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&M.rotation()(0), 9);
}
unsigned int TrajectoryOperationCubic::size() const
{
	return 6;
}

const TrajectorySample & TrajectoryOperationCubic::operator()(double)
{
	return m_sample;
}

const TrajectorySample & TrajectoryOperationCubic::computeNext()
{

	Vector3d rot_diff_vec;
	Vector3d cubic_tra;
	Vector3d cubic_rot_tra;
	Matrix3d cubic_rot;
	m_sample.resize(12, 6);

	Matrix3d rot_diff = m_init.linear().inverse() * m_goal.linear() ;
//	Quaterniond a1(rot_diff);
//	rot_diff_vec = a1.toRotationMatrix().eulerAngles(2, 1, 0);
	//cout << rot_diff_vec.transpose() << endl;


	//for (int i = 0;i < 3;i++) {
	//	//if (abs(rot_diff_vec(i)) > 1.570) {
	//	//	if (rot_diff_vec(i) > 0)
	//	//	rot_diff_vec(i) = 3.1415926535897932384 - rot_diff_vec(i);
	//	//	else
	//	//	rot_diff_vec(i) = rot_diff_vec(i) + 3.1415926535897932384;
	//	//}
	//	if (abs((abs(rot_diff_vec(i)) - 3.1415)) < 0.001)
	//	{
	//		rot_diff_vec(i) = 0.0;
	//	}

	//}
	//cout << rot_diff_vec.transpose() << endl;

	typedef Eigen::Matrix<double, 9, 1> Vector9;
	if (m_time <= m_stime) {
		m_sample.pos.head<3>() = m_init.translation();
		m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&m_init.rotation()(0), 9);

		return m_sample;
	}
	else if (m_time >= m_stime + m_duration) {
		m_sample.pos.head<3>() = m_goal.translation();
		m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&m_goal.rotation()(0), 9);
		return m_sample;
	}
	else {
		double a0, a1, a2, a3;
		double r0, r1, r2, r3;
		for (int i = 0; i < 3; i++) {
			a0 = m_init.translation()(i);
			a1 = 0.0; //m_init.vel(i);
			a2 = 3.0 / pow(m_duration, 2) * (m_goal.translation()(i) - m_init.translation()(i));
			a3 = -1.0 * 2.0 / pow(m_duration, 3) * (m_goal.translation()(i) - m_init.translation()(i));

			cubic_tra(i) = a0 + a1 * (m_time - m_stime) + a2 * pow(m_time - m_stime, 2) + a3 * pow(m_time - m_stime, 3);
		}

		//for (int i = 0; i < 3; i++)
		//{
		//	r0 = 0.0;
		//	r1 = 0.0; //m_init.vel(i);
		//	r2 = 3.0 / pow(m_duration, 2) * (rot_diff_vec(i));
		//	r3 = -1.0 * 2.0 / pow(m_duration, 3) * (rot_diff_vec(i));

		//	cubic_rot_tra(i) = r0 + r1 * (m_time - m_stime) + r2 * pow(m_time - m_stime, 2) + r3 * pow(m_time - m_stime, 3);
		//}
		AngleAxisd A(rot_diff);
		double angle;
		r0 = 0.0;
		r1 = 0.0; //m_init.vel(i);
		r2 = 3.0 / pow(m_duration, 2) * (A.angle());
		r3 = -1.0 * 2.0 / pow(m_duration, 3) * (A.angle());

		angle = r0 + r1 * (m_time - m_stime) + r2 * pow(m_time - m_stime, 2) + r3 * pow(m_time - m_stime, 3);
		//AngleAngle_to_Rot(A.axis(), angle);
		//Matrix3d a;
		//a = AngleAngle_to_Rot(A.axis(), angle);

		cubic_rot = m_init.linear()*AngleAngle_to_Rot(A.axis(), angle);
		m_cubic.translation() = cubic_tra;
		m_cubic.linear() = cubic_rot;//cubic_rot;//m_goal.linear();//
									// 	cout << cubic_rot << endl;
									 //std::cout << GetPhi(cubic_rot, m_goal.linear()) << std::endl;
		//Matrix3d b;
		//m_init.linear()*AngleAxisd(rot_diff_vec(0), Eigen::Vector3d::UnitZ())*AngleAxisd(rot_diff_vec(1), Eigen::Vector3d::UnitY())*AngleAxisd(rot_diff_vec(2), Eigen::Vector3d::UnitX());
		//cout << b << endl;
		//Matrix3d c;
		//c = m_init.linear()* Rotate_with_Z(rot_diff_vec(0))*Rotate_with_Y(rot_diff_vec(1))*Rotate_with_X(rot_diff_vec(2));
		//cout << "goal"<< c << endl;
		m_sample.pos.head<3>() = m_cubic.translation();
		m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&m_cubic.rotation()(0), 9);


		return m_sample;

	}



}

void TrajectoryOperationCubic::getLastSample(TrajectorySample & sample) const
{
	sample = m_sample;
}

bool TrajectoryOperationCubic::has_trajectory_ended() const
{
	return true;
}

void TrajectoryOperationCubic::setGoalSample(Transform3d goal_M)
{
	m_goal = goal_M;
}
void TrajectoryOperationCubic::setInitSample(Transform3d init_M)
{
	m_init = init_M;
}
void TrajectoryOperationCubic::setDuration(const double & duration)
{
	m_duration = duration;
}
void TrajectoryOperationCubic::setCurrentTime(const double & time)
{
	m_time = time;
}
void TrajectoryOperationCubic::setStartTime(const double & time)
{
	m_stime = time;
}

void TrajectoryOperationCubic::setReference(const Transform3d ref) {
	m_sample.resize(12, 6);
	m_sample.pos.head<3>() = ref.translation();
	typedef Eigen::Matrix<double, 9, 1> Vector9;
	m_sample.pos.tail<9>() = Eigen::Map<const Vector9>(&ref.rotation()(0), 9);
}

