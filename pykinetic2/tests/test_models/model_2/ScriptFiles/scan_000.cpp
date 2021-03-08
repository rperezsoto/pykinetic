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
const double k00 = 7.2645705800e+09;
const double k01 = 2.1243982518e+11;
const double k02 = 7.2645705800e+09;
const double k03 = 2.1243982518e+11;
const double k04 = 1.3583198435e-02;
const double k05 = 4.6448966735e-04;
const double k06 = 3.9721718842e-01;
const double k07 = 4.6448966735e-04;
const double k08 = 2.4841851412e+08;
const double k09 = 4.5937881793e+07;
const double k10 = 2.1243982518e+11;
const double k11 = 4.5937881793e+07;

void model( const state_type &x , state_type &dxdt , const double t)
{
    // Ratelaws
    double r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11;
    r00 = k00*x[0]*x[1];
    r01 = k01*x[3];
    r02 = k02*x[2]*x[3];
    r03 = k03*x[4];
    r04 = k04*x[0]*x[1];
    r05 = k05*x[6];
    r06 = k06*x[3];
    r07 = k07*x[6];
    r08 = k08*x[4];
    r09 = k09*x[5];
    r10 = k10*x[5];
    r11 = k11*x[2]*x[6];
    // MassBalances
    dxdt[0] = -r00+r01-r04+r05;
    dxdt[1] = -r00+r01-r04+r05;
    dxdt[2] = -r02+r03+r10-r11;
    dxdt[3] = +r00-r01-r02+r03-r06+r07;
    dxdt[4] = +r02-r03-r08+r09;
    dxdt[5] = +r08-r09-r10+r11;
    dxdt[6] = +r04-r05+r06-r07+r10-r11;
}

//[ublas_main
int main()
{
    outfile.open("scan_000.data", ofstream::out);
    // time parameters
    double tini,tend,tstep,trep;
    tini = 0.0;
    tend = 1E+02; // Final time, s
    tstep = 1E-12; // Timestep, s
    trep = 1E-1; // Timestep to report, s

    // Convergence Parameters
    double rtol,atol;
	rtol=1E-6;atol=1E-12;

    state_type x(7);
    for (boost::numeric::ublas::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i){
        x[*i] = 0.0;
    }
    
    // Initial concentrations: 
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 0.1;

    // Run the solver
    typedef dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<state_type> > > stepper_type;
    integrate_const(stepper_type(), model, x, tini, tend, trep, write_cout);
    outfile.close();
    return 0;
}
//]
