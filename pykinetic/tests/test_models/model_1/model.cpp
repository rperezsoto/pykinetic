#include <iostream>
#include <fstream>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector < double > state_type;
typedef boost::numeric::ublas::matrix < double > matrix_type;

ofstream outfile;

string State2string(const state_type &v, double t)
{
    string s = boost::lexical_cast<std::string>(t);
    s += '\t';
    for (boost::numeric::ublas::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
        s += boost::lexical_cast<std::string>(*i) + '\t';
    return s;
}

void write_cout(const state_type &x, const double t)
{
    string MyVector = State2string(x, t);
    outfile << MyVector << '\n';
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

struct model
{
    void operator()(const state_type &x, state_type &dxdt, const double t)
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
};

struct jacobian
{
    void operator()(const state_type &x, matrix_type &Jac, const double t, state_type &dfdt)
    {
        Jac(0,0) = -k00*x[2]-k04*x[3]-k08*x[4]-k13*x[6]-k14*x[7];
        Jac(0,1) = +k12*x[5]+k15*x[6];
        Jac(0,2) = -k00*x[0];
        Jac(0,3) = +k01-k04*x[0];
        Jac(0,4) = -k08*x[0];
        Jac(0,5) = +k05+k12*x[1];
        Jac(0,6) = +k09-k13*x[0]+k15*x[1];
        Jac(0,7) = -k14*x[0];
        Jac(1,0) = +k13*x[6]+k14*x[7];
        Jac(1,1) = -k02*x[2]-k06*x[3]-k10*x[4]-k12*x[5]-k15*x[6];
        Jac(1,2) = -k02*x[1];
        Jac(1,3) = -k06*x[1];
        Jac(1,4) = +k03-k10*x[1];
        Jac(1,5) = -k12*x[1];
        Jac(1,6) = +k07+k13*x[0]-k15*x[1];
        Jac(1,7) = +k11+k14*x[0];
        Jac(2,0) = -k00*x[2];
        Jac(2,1) = -k02*x[2];
        Jac(2,2) = -k00*x[0]-k02*x[1];
        Jac(2,3) = +k01;
        Jac(2,4) = +k03;
        Jac(2,5) = 0;
        Jac(2,6) = 0;
        Jac(2,7) = 0;
        Jac(3,0) = +k00*x[2]-k04*x[3];
        Jac(3,1) = -k06*x[3];
        Jac(3,2) = +k00*x[0];
        Jac(3,3) = -k01-k04*x[0]-k06*x[1];
        Jac(3,4) = 0;
        Jac(3,5) = +k05;
        Jac(3,6) = +k07;
        Jac(3,7) = 0;
        Jac(4,0) = -k08*x[4];
        Jac(4,1) = +k02*x[2]-k10*x[4];
        Jac(4,2) = +k02*x[1];
        Jac(4,3) = 0;
        Jac(4,4) = -k03-k08*x[0]-k10*x[1];
        Jac(4,5) = 0;
        Jac(4,6) = +k09;
        Jac(4,7) = +k11;
        Jac(5,0) = +k04*x[3]+k13*x[6];
        Jac(5,1) = -k12*x[5];
        Jac(5,2) = 0;
        Jac(5,3) = +k04*x[0];
        Jac(5,4) = 0;
        Jac(5,5) = -k05-k12*x[1];
        Jac(5,6) = +k13*x[0];
        Jac(5,7) = 0;
        Jac(6,0) = +k08*x[4]-k13*x[6]+k14*x[7];
        Jac(6,1) = +k06*x[3]+k12*x[5]-k15*x[6];
        Jac(6,2) = 0;
        Jac(6,3) = +k06*x[1];
        Jac(6,4) = +k08*x[0];
        Jac(6,5) = +k12*x[1];
        Jac(6,6) = -k07-k09-k13*x[0]-k15*x[1];
        Jac(6,7) = +k14*x[0];
        Jac(7,0) = -k14*x[7];
        Jac(7,1) = +k10*x[4]+k15*x[6];
        Jac(7,2) = 0;
        Jac(7,3) = 0;
        Jac(7,4) = +k10*x[1];
        Jac(7,5) = 0;
        Jac(7,6) = +k15*x[1];
        Jac(7,7) = -k11-k14*x[0];
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
        dfdt[5] = 0;
        dfdt[6] = 0;
        dfdt[7] = 0;
    }
};

int main()
{
    outfile.open("model.data", ofstream::out);
    // time parameters
    double tini,tend,trep;
    tini = 0.0;
    tend = 1E-9; // Final time, s
    trep = 1E-11; // s

    // Convergence Parameters
    double rtol,atol;
    rtol=1E-6;atol=1E-12;

    state_type x(8);
    for (boost::numeric::ublas::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i){
        x[*i] = 0.0;
    }

    // Initial concentrations
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 1;

    // Run the solver
    typedef rosenbrock4 < double > stepper_type;
    integrate_const(make_controlled(atol, rtol, stepper_type()),
                    make_pair(model(), jacobian()),
                    x, tini, tend, trep,
                    write_cout);
    outfile.close();
    return 0;
}
