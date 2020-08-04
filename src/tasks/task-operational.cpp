#include "tasks/task-operational.h"
#include "utils/utils.h"

using namespace std;
namespace HQP
{
  namespace tasks
  {
    using namespace constraint;
    using namespace trajectories;

	TaskOperationalSpace::TaskOperationalSpace(const std::string & name, RobotModel & robot, const Index & frameid):
      TaskMotion(name, robot),
      m_frame_id(frameid),
      m_constraint(name, 3, robot.nv()),
      m_ref(12, 6)
    {
      m_v_ref.setZero();
      m_a_ref.setZero();
      m_M_ref.setIdentity();
      m_wMl.setIdentity();
      m_p_error_vec.setZero(6);
      m_v_error_vec.setZero(6);
      m_p.resize(12);
      m_v.resize(6);
      m_p_ref.resize(12);
      m_v_ref_vec.resize(6);
      m_Kp.setZero(6);
      m_Kd.setZero(6);
      m_a_des.setZero(6);
      m_J.setZero(6, robot.nv());
	  m_singular = false;
	  m_beta = 1;
	  transition_flag = false;
	  m_old_sol.setZero(robot.nv());
    }

    int TaskOperationalSpace::dim() const
    {
      return 6;
    }

    const VectorXd & TaskOperationalSpace::Kp() const { return m_Kp; }

    const VectorXd & TaskOperationalSpace::Kd() const { return m_Kd; }

    void TaskOperationalSpace::Kp(Cref_vectorXd Kp)
    {
      assert(Kp.size()==6);
      m_Kp = Kp;
    }

    void TaskOperationalSpace::Kd(Cref_vectorXd Kd)
    {
      assert(Kd.size()==6);
      m_Kd = Kd;
    }
		
    void TaskOperationalSpace::setReference(TrajectorySample & ref)
    {
      m_ref = ref;
      m_M_ref.translation() = ref.pos.head<3>();
      typedef Eigen::Matrix<double, 3, 3> Matrix3;
      m_M_ref.linear()= Eigen::Map<const Matrix3>(&ref.pos(3), 3, 3);

      m_v_ref = MotionVector<double>(ref.vel);
      m_a_ref = MotionVector<double>(ref.acc);
    }

    const TrajectorySample & TaskOperationalSpace::getReference() const
    {
      return m_ref;
    }

    const VectorXd & TaskOperationalSpace::position_error() const
    {
      return m_p_error_vec;
    }

    const VectorXd & TaskOperationalSpace::velocity_error() const
    {
      return m_v_error_vec;
    }

    const VectorXd & TaskOperationalSpace::position() const
    {
      return m_p;
    }

    const VectorXd & TaskOperationalSpace::velocity() const
    {
      return m_v;
    }

    const VectorXd & TaskOperationalSpace::position_ref() const
    {
      return m_p_ref;
    }

    const VectorXd & TaskOperationalSpace::velocity_ref() const
    {
      return m_v_ref_vec;
    }

    const VectorXd & TaskOperationalSpace::getDesiredAcceleration() const
    {
      return m_a_des;
    }

	VectorXd TaskOperationalSpace::getAcceleration(Cref_vectorXd dv) const
    {
      return m_constraint.matrix()*dv + m_drift.vector();
    }

    //Index TaskOperationalSpace::frame_id() const
    //{
    //  return m_frame_id;
    //}

    const ConstraintBase & TaskOperationalSpace::getConstraint() const
    {
      return m_constraint;
    }

    const ConstraintBase & TaskOperationalSpace::compute(const double t, Cref_vectorXd q, Cref_vectorXd v)
    {
         Transform3d oMi;
      MotionVector<double> v_frame;
      
      m_robot.getUpdateKinematics(q, v);
      oMi = m_robot.getTransformation(m_frame_id);

      v_frame = m_robot.getPointVelocity(m_frame_id);
      m_drift.setZero(); // check acc
      m_wMl.linear() = oMi.linear();

      Transform3d b;
      b = m_M_ref.inverse() * oMi;
      m_p_error = log6(b);
      
      m_v_error = v_frame - m_v_ref; //check actInv
      m_v_error = actinv(oMi, m_v_error.vector());
      m_p_error_vec = m_p_error.vector();
      m_v_error_vec = m_v_error.vector();

      m_p_ref.head<3>() = m_M_ref.translation();
      typedef Eigen::Matrix<double, 9, 1> Vector9;
      m_p_ref.tail<9>() = Eigen::Map<const Vector9>(&m_M_ref.rotation()(0), 9);
      
      m_v_ref_vec = m_v_ref.vector();
      
      m_p.head<3>() = oMi.translation();
      m_p.tail<9>() = Eigen::Map<const Vector9>(&oMi.rotation()(0), 9);

      m_v = v_frame.vector();

      m_a_des = -m_Kp.cwiseProduct(m_p_error.vector())
            - m_Kd.cwiseProduct(m_v_error.vector());
      //+ m_a_ref.actInv(m_wMl).vector();

      Vector3d delphi = GetPhi(oMi.linear(), m_M_ref.linear());
      Vector6d m_a_des_new1;
      m_a_des_new1.head(3) = 1000 * (m_M_ref.translation() - oMi.translation()) - 2.0 * sqrt(600.0) * v_frame.linear();
      m_a_des_new1.tail(3) = -200.0 * delphi - 30.0*v_frame.angular();

      m_J = m_robot.getJacobian(m_frame_id);

      m_constraint.setMatrix(m_J);
      if (!transition_flag)
        m_constraint.setVector(m_beta * m_a_des_new1);
      else
        m_constraint.setVector(m_beta * m_a_des_new1 + (1 - m_beta) * m_J * m_old_sol);
      return m_constraint;

      }   
  }
}
