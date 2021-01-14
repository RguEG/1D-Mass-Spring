//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
	Eigen::VectorXd q_1,q_2,q_3,q_4;
	Eigen::VectorXd q_dot_1, q_dot_2, q_dot_3, q_dot_4;
	Eigen::VectorXd xforce;
	force(xforce, q, qdot);
	q_dot_1 = xforce / mass;
	q_dot_2 = dt / 2 * q_dot_1;
	force(xforce, q, q_dot_2);
	q_dot_2 = xforce / mass;
	q_dot_3 = dt / 2 * q_dot_2;
	force(xforce, q, q_dot_3);
	q_dot_3 = xforce / mass;
	q_dot_4 = dt / 2 * q_dot_3;
	force(xforce, q, q_dot_4);
	q_dot_4 = xforce / mass;
	q_dot = q_dot + dt / 6 * (q_dot_1 + 2 * q_dot_2 + 2 * q_dot_3 + q_dot_4);
	q = q + dt * qdot;
}