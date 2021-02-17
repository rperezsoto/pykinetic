#include <iostream>
#include <fstream>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double > state_type;

ofstream outfile;

string State2string(const state_type &v ,double t)
{
    string s = boost::lexical_cast<std::string>(t);
    s += '\t';
    for (boost::numeric::ublas::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
        s += boost::lexical_cast<std::string>(*i) + '\t';
    return s;
}

void write_cout(const state_type &x, const double t)
{
    string MyVector = State2string(x,t);
    outfile <<  MyVector << '\n';
}

// Kinetic Constants
const double k00 = 2.1243982516e+11;
const double k01 = 3.2362038820e-39;
const double k02 = 2.1243982516e+11;
const double k03 = 3.2362038820e-39;
const double k04 = 2.1243982521e+11;
const double k05 = 2.6769180696e-24;
const double k06 = 2.1243982521e+11;
const double k07 = 6.6268592390e-26;
const double k08 = 2.1243982521e+11;
const double k09 = 6.6268592390e-26;
const double k10 = 2.1243982521e+11;
const double k11 = 2.6769180696e-24;
const double k12 = 2.1243982521e+11;
const double k13 = 5.2590657658e+09;
const double k14 = 2.1243982521e+11;
const double k15 = 5.2590657658e+09;

void model( const state_type &x , state_type &dxdt , const double t)
{
    // Ratelaws
    double r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11,r12,r13,r14,r15;
    r00 = k00*x[0]*x[2];
    r01 = k01*x[3];
    r02 = k02*x[1]*x[2];
    r03 = k03*x[4];
    r04 = k04*x[0]*x[3];
    r05 = k05*x[5];
    r06 = k06*x[1]*x[3];
    r07 = k07*x[6];
    r08 = k08*x[0]*x[4];
    r09 = k09*x[6];
    r10 = k10*x[1]*x[4];
    r11 = k11*x[7];
    r12 = k12*x[1]*x[5];
    r13 = k13*x[0]*x[6];
    r14 = k14*x[0]*x[7];
    r15 = k15*x[1]*x[6];
    // MassBalances
    dxdt[0] = -r00+r01-r04+r05-r08+r09+r12-r13-r14+r15;
    dxdt[1] = -r02+r03-r06+r07-r10+r11-r12+r13+r14-r15;
    dxdt[2] = -r00+r01-r02+r03;
    dxdt[3] = +r00-r01-r04+r05-r06+r07;
    dxdt[4] = +r02-r03-r08+r09-r10+r11;
    dxdt[5] = +r04-r05-r12+r13;
    dxdt[6] = +r06-r07+r08-r09+r12-r13+r14-r15;
    dxdt[7] = +r10-r11-r14+r15;
}

//[ublas_main
int main()
{
    outfile.open("Main.data", ofstream::out);
	// time parameters
    double tini,tend,tstep,trep;
    tini = 0.0;
    tend = 1E-9; // Final time, s
    tstep = 1E-12; // Timestep, s
	trep = 1E-11; // Timestep to report, s

	// Convergence Parameters
	double rtol,atol;
	rtol=1E-6;atol=1E-12;

	state_type x(8);
	for (boost::numeric::ublas::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i){
        x[*i] = 0.0;
    }
	x[0] = 0.5;
	x[1] = 0.5;
	x[2] = 1;
    typedef dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<state_type> > > stepper_type;
    integrate_const(stepper_type(), model, x, tini, tend, trep, write_cout);
    outfile.close();
    return 0;
}
//]
