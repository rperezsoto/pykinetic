#include <iostream>
#include <fstream>


#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double > state_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


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

struct model
{
    void operator()( const state_type &x , state_type &dxdt , const double t)
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
};

struct jacobian
{
    void operator()(const state_type &x, matrix_type &Jac, const double t, state_type &dfdt)
    {
        // Non-zero Elements of the Jacobian
        Jac(0,0) = -k00*x[1]-k04*x[1];
        Jac(0,1) = -k00*x[0]-k04*x[0];
        Jac(0,2) = 0;
        Jac(0,3) = +k01;
        Jac(0,4) = 0;
        Jac(0,5) = 0;
        Jac(0,6) = +k05;
        Jac(1,0) = -k00*x[1]-k04*x[1];
        Jac(1,1) = -k00*x[0]-k04*x[0];
        Jac(1,2) = 0;
        Jac(1,3) = +k01;
        Jac(1,4) = 0;
        Jac(1,5) = 0;
        Jac(1,6) = +k05;
        Jac(2,0) = 0;
        Jac(2,1) = 0;
        Jac(2,2) = -k02*x[3]-k11*x[6];
        Jac(2,3) = -k02*x[2];
        Jac(2,4) = +k03;
        Jac(2,5) = +k10;
        Jac(2,6) = -k11*x[2];
        Jac(3,0) = +k00*x[1];
        Jac(3,1) = +k00*x[0];
        Jac(3,2) = -k02*x[3];
        Jac(3,3) = -k01-k02*x[2]-k06;
        Jac(3,4) = +k03;
        Jac(3,5) = 0;
        Jac(3,6) = +k07;
        Jac(4,0) = 0;
        Jac(4,1) = 0;
        Jac(4,2) = +k02*x[3];
        Jac(4,3) = +k02*x[2];
        Jac(4,4) = -k03-k08;
        Jac(4,5) = +k09;
        Jac(4,6) = 0;
        Jac(5,0) = 0;
        Jac(5,1) = 0;
        Jac(5,2) = +k11*x[6];
        Jac(5,3) = 0;
        Jac(5,4) = +k08;
        Jac(5,5) = -k09-k10;
        Jac(5,6) = +k11*x[2];
        Jac(6,0) = +k04*x[1];
        Jac(6,1) = +k04*x[0];
        Jac(6,2) = -k11*x[6];
        Jac(6,3) = +k06;
        Jac(6,4) = 0;
        Jac(6,5) = +k10;
        Jac(6,6) = -k05-k07-k11*x[2];
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
        dfdt[5] = 0;
        dfdt[6] = 0;
    }    
};

//[ublas_main
int main()
{
    outfile.open("model.data", ofstream::out);
    // time parameters
    double tini,tend,tstep,trep;
    tini = 0.0;
    tend = 1E+02; // Final time, s
    tstep = 1E-12; // Timestep, s
    trep = 1E-1; // s

    // Convergence Parameters
    double rtol,atol;
	rtol=1E-6;atol=1E-12;

    state_type x(7);
    for (boost::numeric::ublas::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i){
        x[*i] = 0.0;
    }

    // Initial concentrations
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 0.1;

    // Run the solver
    typedef rosenbrock4< double > stepper_type;
    integrate_const(make_controlled(atol, rtol, stepper_type()),
                    make_pair(model(), jacobian()),
                    x, tini, tend, trep, 
                    write_cout);
    outfile.close();
    return 0;
}
//]
