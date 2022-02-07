#include<bits/stdc++.h>

using namespace std;

class useful_functions{ //this class contains the christoffel coefficients and the length scales Delta, Sigma, and A as defined https://en.wikipedia.org/wiki/Kerr_metric#Metric 
    
public:
    
    //first the length scales

    double Kerr_Delta(double r,double a){
        double kerr_delta=r*r-2*r+a*a;
        return kerr_delta;
    }
    
    double Kerr_Sigma(double r,double a,double theta){
        double kerr_sigma=r*r+a*a*cos(theta)*cos(theta);
        return kerr_sigma;
    }
    
    double Kerr_A(double r,double a,double theta){
        double kerr_a=(r*r+a*a)*(r*r+a*a)-a*a*Kerr_Delta(r,a)*sin(theta)*sin(theta);
        return kerr_a;
    }

    //These are the coefficients of the line element in Boyer-Lindquist coordinates
    
    double line_element_coeff_1(double r,double a,double theta){
        double LINE_ELEMENT_COEFF_1=-(1-2*r/Kerr_Sigma(r,a,theta));
        return LINE_ELEMENT_COEFF_1;
    }
    
    double line_element_coeff_2(double r,double a,double theta){
        double LINE_ELEMENT_COEFF_2=-2*a*r*sin(theta)*sin(theta)/Kerr_Sigma(r,a,theta);
        return LINE_ELEMENT_COEFF_2;
    }
    
    double line_element_coeff_3(double r,double a,double theta){
        double LINE_ELEMENT_COEFF_3=Kerr_Sigma(r,a,theta)/Kerr_Delta(r,a);
        return LINE_ELEMENT_COEFF_3;
    }
    
    double line_element_coeff_4(double r,double a,double theta){
        double LINE_ELEMENT_COEFF_4=Kerr_Sigma(r,a,theta);
        return LINE_ELEMENT_COEFF_4;
    }
    
    double line_element_coeff_5(double r,double a,double theta){
        double LINE_ELEMENT_COEFF_5=(r*r+a*a+ 2*a*a*r*sin(theta)*sin(theta)/Kerr_Sigma(r,a,theta))*sin(theta)*sin(theta);
        return LINE_ELEMENT_COEFF_5;
    }
    
    //now the christoffel coefficients
    
    //the radial coefficients are first, the angular ones second
    
    
    double gamma_r_tt(double r,double a,double theta){
        double GAMMA_R_TT=Kerr_Delta(r,a)*(r*r-a*a*cos(theta)*cos(theta))/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_R_TT;
    }
    
    double gamma_r_rr(double r,double a,double theta){
        double GAMMA_R_RR=(r*a*a*sin(theta)*sin(theta)-(r*r-a*a*cos(theta)*cos(theta)))/(Kerr_Sigma(r,a,theta)*Kerr_Delta(r,a));
        return GAMMA_R_RR;
    }
    
    double gamma_r_oo(double r,double a,double theta){
        double GAMMA_R_OO=- r*Kerr_Delta(r,a)/Kerr_Sigma(r,a,theta);
        return GAMMA_R_OO;
    }
    
    double gamma_r_pp(double r,double a,double theta){
        double GAMMA_R_PP= Kerr_Delta(r,a)*sin(theta)*sin(theta)*(-r*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)+a*a*sin(theta)*sin(theta)*(r*r-a*a*cos(theta)*cos(theta)))/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_R_PP;
    }
    
    double gamma_r_tp(double r,double a,double theta){
        double GAMMA_R_TP=- Kerr_Delta(r,a)*a*sin(theta)*sin(theta)*(r*r-a*a*cos(theta)*cos(theta))/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_R_TP;
    }
    
    double  gamma_r_ro(double r,double a,double theta){
        double GAMMA_R_RO=- a*a*sin(theta)*cos(theta)/(Kerr_Sigma(r,a,theta));
        return GAMMA_R_RO;
    }
    
    double  gamma_o_tt(double r,double a,double theta){
        double GAMMA_O_TT=-2*a*a*r*sin(theta)*cos(theta)/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_O_TT;
    }
    
    double gamma_o_rr(double r,double a,double theta){
        double GAMMA_O_RR=a*a*sin(theta)*cos(theta)/(Kerr_Sigma(r,a,theta)*Kerr_Delta(r,a));
        return GAMMA_O_RR;
    }
    
    double gamma_o_oo(double r,double a,double theta){
        double GAMMA_O_OO=- a*a*sin(theta)*cos(theta)/Kerr_Sigma(r,a,theta);
        return GAMMA_O_OO;
    }
    
    double gamma_o_pp(double r,double a,double theta){
        double GAMMA_O_PP=- sin(theta)*cos(theta)*(Kerr_A(r,a,theta)*Kerr_Sigma(r,a,theta)+2*(r*r+a*a)*a*a*r*sin(theta)*sin(theta))/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_O_PP;
    }
    
    double gamma_o_tp(double r,double a,double theta){
        double GAMMA_O_TP=2*a*r*(r*r+a*a)*sin(theta)*cos(theta)/(Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta)*Kerr_Sigma(r,a,theta));
        return GAMMA_O_TP;
    }
    
    double gamma_o_ro(double r,double a,double theta){
        double GAMMA_O_RO=r/Kerr_Sigma(r,a,theta);
        return GAMMA_O_RO;
    }
    
    //to catch errors propagating through the integration before they do too much damage, this function is defined. The flag should always be -1 for a photon!
    
    double check_integral_error_propagation(double r,double a,double theta,double td,double rd,double od,double pd){
        double propagated_integral_error_flag=(line_element_coeff_3(r,a,theta)*rd*rd+line_element_coeff_5(r,a,theta)*pd*pd+line_element_coeff_4(r,a,theta)*od*od+2*line_element_coeff_2(r,a,theta)*td*pd)/(line_element_coeff_1(r,a,theta)*td*td);
        return propagated_integral_error_flag;
    }
    
    //there are two quantities conserved along the geodesic - the total photon energy and the angular momentum. (if a photon hits the black hole the angular momentum gets transferred!)
    
    double conserved_energy(double r,double a,double theta,double td,double pd){
        double CONSERVED_ENERGY=- line_element_coeff_1(r,a,theta)*td-line_element_coeff_2(r,a,theta)*pd;
        return CONSERVED_ENERGY;
    }
    
    double conserved_angular_momentum(double r,double a,double theta,double td,double pd){
        double CONSERVED_ANGULAR_MOMENTUM=line_element_coeff_5(r,a,theta)*pd+line_element_coeff_2(r,a,theta)*td ;
        return CONSERVED_ANGULAR_MOMENTUM;
    }
    
    double normalized_angular_momentum(double r,double a,double theta,double td,double pd){
        double NORMALIZED_ANGULAR_MOMENTUM=conserved_angular_momentum(r,a,theta,td,pd)/conserved_energy(r,a,theta,td,pd);
        return NORMALIZED_ANGULAR_MOMENTUM;
    }
};



class Initial_conditions{
    
public:
    
    //these functions take a position in the image plane - say initial_x,initial_y - and have an angle initial_theta from the center of the black hole. 
    //they then turn them into initial conditions in Boyer-Lindquist coordinates
    //also the initial photon 4-velocity is generated

    //first initialize the position
    
    useful_functions uf;
    
    double initial_r(double d,double initial_x,double initial_y){
        double R_INIT=sqrt(d*d+initial_x*initial_x+initial_y*initial_y);
        return R_INIT;
    }
    
    double initial_theta(double d,double initial_x,double initial_y,double initial_theta){
        double THETA_INIT=acos( (d*cos(initial_theta)+ initial_y*sin(initial_theta))/initial_r(d,initial_x,initial_y) );
        return THETA_INIT;
    }
    
    double initial_phi(double d,double initial_x,double initial_y,double initial_theta){
        double PHI_INIT=atan(initial_x/(d*sin(initial_theta)-initial_y*cos(initial_theta)));
        return PHI_INIT;
    }
    
    // now initialize the 4-velocity
    
    double initial_vr(double d,double initial_x,double initial_y){
        double VR_INIT=-d/initial_r(d,initial_x,initial_y);
        return VR_INIT;
    }
    
    double initial_vo(double d,double initial_x,double initial_y,double initial_theta){
        double VO_INIT=-(-cos(initial_theta)+d/(pow(initial_r(d,initial_x,initial_y),2))*(d*cos(initial_theta)+initial_y*sin(initial_theta)))/sqrt(pow(initial_r(d,initial_x,initial_y),2)-pow((d*cos(initial_theta)+initial_y*sin(initial_theta)),2) );
        return VO_INIT;
    }
    
    double initial_vp(double d,double initial_x,double initial_y,double initial_theta){
        double VP_INIT=- (-initial_x*sin(initial_theta)/(pow(d*sin(initial_theta)-initial_y*cos(initial_theta),2)+initial_x*initial_x));
        return VP_INIT;
    }
    
    double initial_vt(double a,double d,double initial_x,double initial_y,double initial_theta){
        double creative_variable_1=uf.line_element_coeff_1(initial_r(d,initial_x,initial_y),a,initial_theta);
        double creative_variable_2=2*uf.line_element_coeff_2(initial_r(d,initial_x,initial_y),a,initial_theta)*initial_vp(d,initial_x,initial_y,initial_theta);
        double creative_variable_3 =uf.line_element_coeff_3(initial_r(d,initial_x,initial_y),a,initial_theta)*initial_vr(d,initial_x,initial_y)*initial_vr(d,initial_x,initial_y)+uf.line_element_coeff_4(initial_r(d,initial_x,initial_y),a,initial_theta)*initial_vo(d,initial_x,initial_y,initial_theta)*initial_vo(d,initial_x,initial_y,initial_theta)+uf.line_element_coeff_5(initial_r(d,initial_x,initial_y),a,initial_theta)*initial_vp(d,initial_x,initial_y,initial_theta)*initial_vp(d,initial_x,initial_y,initial_theta);
        
        //for an actual physical verobjectsion, this has to be greater than 0
        double VT_INIT=(-creative_variable_2-sqrt(creative_variable_2*creative_variable_2-4*creative_variable_1*creative_variable_3))/(2*creative_variable_1);
        return VT_INIT;
    }
};

class motion_equations{
    
    
public:
    
    //now it's time for the equations of motion
    //phi and t evolve according to first-order equations
    //r and theta evolve according to second-order geodesic equations
    
    //for the second-order equations we solve them as first order to get velocities and then use those velocities in further calculations
    
    useful_functions uf;
    
    double t_dot(double r,double a,double theta,double j){
        double T_DOT=-(uf.line_element_coeff_5(r,a,theta)+j*uf.line_element_coeff_2(r,a,theta))/(uf.line_element_coeff_5(r,a,theta)*uf.line_element_coeff_1(r,a,theta)-uf.line_element_coeff_2(r,a,theta)*uf.line_element_coeff_2(r,a,theta));
        return T_DOT;
    }
    
    double p_dot(double r,double a,double theta,double j){
        double P_DOT=(j*uf.line_element_coeff_1(r,a,theta)+uf.line_element_coeff_2(r,a,theta))/(uf.line_element_coeff_5(r,a,theta)*uf.line_element_coeff_1(r,a,theta)-uf.line_element_coeff_2(r,a,theta)*uf.line_element_coeff_2(r,a,theta));
        return P_DOT;
    }
    
    double vr_dot(double r,double a,double theta,double td,double rd,double od,double pd){
        double VR_DOT=(-uf.gamma_r_tt(r,a,theta)*td*td-uf.gamma_r_rr(r,a,theta)*rd*rd-uf.gamma_r_oo(r,a,theta)*od*od-uf.gamma_r_pp(r,a,theta)*pd*pd-2*uf.gamma_r_tp(r,a,theta)*td*pd-2*uf.gamma_r_ro(r,a,theta)*rd*od);
        return VR_DOT;
    }
    
    double vo_dot(double r,double a,double theta,double td,double rd,double od,double pd){
        double VO_DOT=-uf.gamma_o_tt(r,a,theta)*td*td-uf.gamma_o_rr(r,a,theta)*rd*rd-uf.gamma_o_oo(r,a,theta)*od*od-uf.gamma_o_pp(r,a,theta)*pd*pd-2*uf.gamma_o_tp(r,a,theta)*td*pd-2*uf.gamma_o_ro(r,a,theta)*rd*od;
        return VO_DOT;
    }
    
};


class Integrators{
    
    //DO NOT TOUCH THIS CODE - MESSING WITH THIS WILL BREAK THE INTEGRATOR IN UNRECOVERABLE WAYS
    //if you want funny results then change the denominators of the RK4 solver
    
public:

    motion_equations motion_equation;
    useful_functions uf;
    
    void rk4_1(double r,double a,double theta,double td,double rd,double od,double pd,double j,double h,double* k1){
        k1[0]=motion_equation.vr_dot(r,a,theta,td,rd,od,pd);
        k1[1]=motion_equation.vo_dot(r,a,theta,td,rd,od,pd);
        k1[2]=motion_equation.t_dot(r,a,theta,j);
        k1[3]=rd;
        k1[4]=od;
        k1[5]=motion_equation.p_dot(r,a,theta,j);
    }
    
    void rk4_2(double r,double a,double theta,double td,double rd,double od,double pd,double j,double h,double* k1,double* k2){
        k2[0]=motion_equation.vr_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,motion_equation.t_dot(r+ h*k1[3]/2,a,theta+h*k1[4]/2,j ),rd+h*k1[0]/2,od+h*k1[1]/2,motion_equation.p_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,j ));
        k2[1]=motion_equation.vo_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,motion_equation.t_dot(r+ h*k1[3]/2,a,theta+h*k1[4]/2,j ),rd+h*k1[0]/2,od+h*k1[1]/2,motion_equation.p_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,j ));
        k2[2]=motion_equation.t_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,j );
        k2[3]=rd+h*k1[0]/2;
        k2[4]=od+h*k1[1]/2;
        k2[5]=motion_equation.p_dot(r+h*k1[3]/2,a,theta+h*k1[4]/2,j );
    }
    
    void rk4_3(double r,double a,double theta,double td,double rd,double od,double pd,double j,double h,double* k2,double* k3){
        k3[0]=motion_equation.vr_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,motion_equation.t_dot(r+ h*k2[3]/2,a,theta+h*k2[4]/2,j ),rd+h*k2[0]/2,od+h*k2[1]/2,motion_equation.p_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,j ));
        k3[1]=motion_equation.vo_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,motion_equation.t_dot(r+ h*k2[3]/2,a,theta+h*k2[4]/2,j ),rd+h*k2[0]/2,od+h*k2[1]/2,motion_equation.p_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,j ));
        k3[2]=motion_equation.t_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,j );
        k3[3]=rd+h*k2[0]/2;
        k3[4]=od+h*k2[1]/2;
        k3[5]=motion_equation.p_dot(r+h*k2[3]/2,a,theta+h*k2[4]/2,j );
    }
    
    void rk4_4(double r,double a,double theta,double td,double rd,double od,double pd,double j,double h,double* k3,double* k4){
        k4[0]=motion_equation.vr_dot(r+h*k3[3],a,theta+h*k3[4],motion_equation.t_dot(r+ h*k3[3],a,theta+h*k3[4],j ),rd+h*k3[0],od+h*k3[1],motion_equation.p_dot(r+h*k3[3],a,theta+h*k3[4],j ));
        k4[1]=motion_equation.vo_dot(r+h*k3[3],a,theta+h*k3[4],motion_equation.t_dot(r+ h*k3[3],a,theta+h*k3[4],j ),rd+h*k3[0],od+h*k3[1],motion_equation.p_dot(r+h*k3[3],a,theta+h*k3[4],j ));
        k4[2]=motion_equation.t_dot(r+h*k3[3],a,theta+h*k3[4],j );
        k4[3]=rd+h*k3[0];
        k4[4]=od+h*k3[1];
        k4[5]=motion_equation.p_dot(r+h*k3[3],a,theta+h*k3[4],j );
    }
    
    void one_time_step(double* Y,double a,double j,double h,double* k1,double* k2,double* k3,double* k4){
        double r=Y[3];
        double theta=Y[4];
        double rd=Y[0];
        double od=Y[1];
        
        double td=motion_equation.t_dot(r,a,theta,j);
        double pd=motion_equation.p_dot(r,a,theta,j);
        
        rk4_1(r,a,theta,td,rd,od,pd,j,h,k1);
        rk4_2(r,a,theta,td,rd,od,pd,j,h,k1,k2);
        rk4_3(r,a,theta,td,rd,od,pd,j,h,k2,k3);
        rk4_4(r,a,theta,td,rd,od,pd,j,h,k3,k4);
        
        Y[0]=Y[0]+1.0/6.0*h*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        Y[1]=Y[1]+1.0/6.0*h*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        Y[2]=Y[2]+1.0/6.0*h*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        Y[3]=Y[3]+1.0/6.0*h*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);
        Y[4]=Y[4]+1.0/6.0*h*(k1[4]+2*k2[4]+2*k3[4]+k4[4]);
        Y[5]=Y[5]+1.0/6.0*h*(k1[5]+2*k2[5]+2*k3[5]+k4[5]);
    }
    
    double next_time_step(double initial_x,double rad,double vrad,double thet,double vthet,double ph,double vph,double F){
        // gets the next time-step in the algorithm.
        double elm[3];
        double min;
        double dT;
        if (initial_x != 0){// dT is set by F multiplied by the shortest variability time-scale in the problem.
            elm[0]=abs(rad/vrad);
            elm[1]=abs(thet/vthet);
            elm[2]=abs(ph/vph);
            min=elm[0];
            for (int m=0; m<2; m ++){
                if (elm[m]<min){
                    min=elm[m];
                }
            }
            dT=F*min;
        }
        else{//If initial_x=0,then \phi=0,and we will pass over the pole \theta=0,at some point in the trajectory.
            // We must only consider the radial time-scale in that case.
            dT=F*abs(rad/vrad);
        }
    
        return dT;
    }
};

class photon{
    
  // This will solve the geodesic equations for a single photon
    
public:
    useful_functions uf;
    motion_equations motion_equation;
    Integrators integrator;
    Initial_conditions start;
    
    void raytracing_output(vector<double> arr,int len_arr){
        // Takes a vector input arr,with length len_arr
        // and outputs it as one big string.
        //
        // This string can be converted into a numpy array
        // quite easily within python.
        stringstream ss;
        for(int i=0; i<len_arr; ++i)
        {
          if(i != 0)
            ss<<",";
          ss<<arr[i];
        }
        string s=ss.str();
        cout<<s<<endl;
    }
    
    double innermost_stable_circular_orbit_prograde(double a){
        // Returns the ISCO radius for prograde spins,a>0.
        double Z_1=1+pow(1-a*a,1.0/3.0)*(pow(1+a,1.0/3.0)+pow(1-a,1.0/3.0));
        double Z_2=sqrt(3*a*a+Z_1*Z_1);
        double innermost_stable_circular_orbit= (3+Z_2- sqrt((3-Z_1)*(3+Z_1+2*Z_2)));
        return innermost_stable_circular_orbit;
    }

    double innermost_stable_circular_orbit_retrograde(double a){
        // Returns the ISCO radius for retrograde spins,a<0.
        double Z_1=1+pow(1-a*a,1.0/3.0)*(pow(1+a,1.0/3.0)+pow(1-a,1.0/3.0));
        double Z_2=sqrt(3*a*a+Z_1*Z_1);
        double innermost_stable_circular_orbit= (3+Z_2+sqrt((3-Z_1)*(3+Z_1+2*Z_2)));
        return innermost_stable_circular_orbit;
    }

    double get_innermost_stable_circular_orbit(double a){
    // Gets the ISCO radius for a given BH spin parameter
        double ri=1.0;
        if (a == 0){
            ri=6;// Schwarzschild result.
        }
        else if (a>0){
            ri=innermost_stable_circular_orbit_prograde(a);// Analytic formula for the ISCO depends on the sign of the spin parameter.
        }
        else {
            ri=innermost_stable_circular_orbit_retrograde(a);// Analytic formula for the ISCO depends on the sign of the spin parameter.
        }
        return ri;
    }
    
    double obtain_event_horizon(double a){
        // Gets the event horizon for a given BH spin parameter.
        double rH= 1+sqrt(1-a*a);
        return rH;
    }

    void compute_orbit(double a,double initial_x,double initial_y,double initial_theta,double d,double F,const int n,double& rf,double& phif,double& thetaf,string method){
        // First we initialise the photons 4-velocities.
        
        double vt0=start.initial_vt(a,d,initial_x,initial_y,initial_theta);
        double vr0=start.initial_vr(d,initial_x,initial_y);
        double vo0=start.initial_vo(d,initial_x,initial_y,initial_theta);
        double vp0=start.initial_vp(d,initial_x,initial_y,initial_theta);
        
        // then we initialise the photons position
        
        double t0=1;
        double r0=start.initial_r(d,initial_x,initial_y);
        double theta0=start.initial_theta(d,initial_x,initial_y,initial_theta);
        double phi0=start.initial_phi(d,initial_x,initial_y,initial_theta);
        
        // calculate the photons normalised angular momentum
        // the important quantity for red-shift calculations
        
        double j=uf.normalized_angular_momentum(r0,a,theta0,vt0,vp0);
        
        // calculate initial r and theta acceleration
        
        double vr_dot0=motion_equation.vr_dot(r0,a,theta0,vt0,vr0,vo0,vp0);
        double vo_dot0=motion_equation.vo_dot(r0,a,theta0,vt0,vr0,vo0,vp0);
        
        // important vectors for algorithm implimentation
        
        vector<double> t(1);
        double dT=0;
        
        // Varying time-stepping,time-step in integration is
        // a fixed fraction of the time-scale of the
        // rapidest changing variable
        
        // d tau == F*min( r/r_dot,theta/theta_dot,phi/phi_dot )
        
        dT=integrator.next_time_step(initial_x,r0,vr0,theta0,vo0,phi0,vp0,F);
        t[0]=dT;
        
        // set-up for the algorithm implimentation
        
        double rk4y[6];
        double k1[6];
        double k2[6];
        double k3[6];
        double k4[6];
        
        
        rk4y[0]=vr0;
        rk4y[1]=vo0;
        rk4y[2]=t0;
        rk4y[3]=r0;
        rk4y[4]=theta0;
        rk4y[5]=phi0;
        
        
        k1[0]=k1[1]=k1[2]=k1[3]=k1[4]=k1[5]=0;
        k2[0]=k2[1]=k2[2]=k2[3]=k2[4]=k2[5]=0;
        k3[0]=k3[1]=k3[2]=k3[3]=k3[4]=k3[5]=0;
        k4[0]=k4[1]=k4[2]=k4[3]=k4[4]=k4[5]=0;
        
        // Quantities calculated at every time-step
        
        vector<double> r(1);
        vector<double> theta(1);
        vector<double> phi(1);


        r[0]=r0;
        theta[0]=theta0;
        phi[0]=phi0;

        int m=0;// For watching the time-steps.
        
        // First method,all photon paths terminated when they move through the disk plane (z=0).
        // Even if they pass between the ISCO and event horizon.
        if (method == "Simple"){
            double rH=obtain_event_horizon(a);
            while (rk4y[3]*cos(rk4y[4])>0){ // stops when z<0
                if (m<n){ // unless algorithm has got stuck (generally happens when photon hits black-hole)
                    if (rk4y[3]>rH){ // if r<rH then  in the black-hole
                        integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                        r.push_back(rk4y[3]);
                        theta.push_back(rk4y[4]); // updates the variables
                        phi.push_back(rk4y[5]);
                        dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                        t.push_back(t[m]+dT);
                    }
                    else{//I am in the black hole.
                        break;
                    }
                    ++m; // in case algorithm gets stuck
                }
                else{// Too many steps.
                    break;
                }
            }
        }
        
        // Second method,allows light rays to pass between the disk and event horizon: rH<r<rI.
        if (method == "Disk"){
            double rI=get_innermost_stable_circular_orbit(a);
            double rH=obtain_event_horizon(a);
            while (rk4y[3]*cos(rk4y[4])>-20){
                if (m<n){
                    if  (rk4y[0]*cos(rk4y[4])-rk4y[3]*sin(rk4y[4])*rk4y[1]<0){
                        if (rk4y[3]>rH){
                            if (abs(rk4y[3]*cos(rk4y[4]))>0.05){
                                    integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                                    t.push_back(t[m]+dT);
                                    ++m;
                            }
                            else if (rk4y[3]<rI){
                                    integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                                    t.push_back(t[m]+dT);
                                    ++m;
                                    }
                            else{
                                break;
                            }
                        }
                        else{
                            break;
                        }
                    }
                    else if (rk4y[3]*cos(rk4y[4])<20){
                        if (rk4y[3]>rH){
                            if (abs(rk4y[3]*cos(rk4y[4]))>0.1){
                                    integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                                    t.push_back(t[m]+dT);
                                    ++m;
                            }
                                else if (rk4y[3]<rI){
                                    integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                                    t.push_back(t[m]+dT);
                                    ++m;
                                }
                                else{
                                    break;
                                }
                        }
                    else{
                        break;
                    }
                    }
                else{
                    break;
                }
            

            }
            }
        }
        
        // Third method,assumes no disk in the z=0 plane.
        if (method == "NoDisk"){
            double rH=obtain_event_horizon(a);
            while ( (rk4y[3]*cos(rk4y[4])>-20) && (rk4y[3]*sin(rk4y[4])*cos(rk4y[5])>-30)){// stops when photon reaches z=-20 or x=-30,or x>30 and v_x>0.
                if ( ((rk4y[0]*sin(rk4y[4])*cos(rk4y[5])+rk4y[3]*rk4y[1]*cos(rk4y[4])*cos(rk4y[5])-rk4y[3]*motion_equation.p_dot(rk4y[3],a,rk4y[4],j)*sin(rk4y[4])*sin(rk4y[5])>0)  && rk4y[3]*sin(rk4y[4])*cos(rk4y[5])>30) ){
                    break;// Breaks if x>30 and v_x>0.
                }
                if (m<n){
                    if  (rk4y[0]*cos(rk4y[4])-rk4y[3]*sin(rk4y[4])*rk4y[1]<0){// if z-velocity is negative.
                        if (rk4y[3]>rH){
                            integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                            r.push_back(rk4y[3]);
                            theta.push_back(rk4y[4]); // updates the variables
                            phi.push_back(rk4y[5]);
                            dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                            t.push_back(t[m]+dT);
                            ++m;
                        }
                        else{
                            break;
                        }
                    }
                    else if (rk4y[3]*cos(rk4y[4])<20){// else I am going up,and z <+20.
                        if (rk4y[3]>rH){
                            integrator.one_time_step(rk4y,a,j,dT,k1,k2,k3,k4); // implements the RK4 algorithm
                            r.push_back(rk4y[3]);
                            theta.push_back(rk4y[4]); // updates the variables
                            phi.push_back(rk4y[5]);
                            dT=integrator.next_time_step(initial_x,rk4y[3],rk4y[0],rk4y[4],rk4y[1],rk4y[5],motion_equation.p_dot(rk4y[3],a,rk4y[4],j),F);
                            t.push_back(t[m]+dT);
                            ++m;
                        }
                    else{
                        break;
                    }
                    }
                else{
                    break;
                }
            

            }
            }
        }

        
        raytracing_output(r,m);
        raytracing_output(theta,m);
        raytracing_output(phi,m);
        raytracing_output(t,m);
        
        cout<<j<<endl;
        
        rf=r[m];
        thetaf=theta[m]; // Edits these variables in the main script
        phif=phi[m];
        
        
    }
};

int main(int argc,const char*argv[]) {
    // We import the class of one photon orbit integrators
    
    photon example_photon;
    
    if (argc != 6){// Should get 1: C++ file name (always happens),2-5: Parameters alphinitial_x,betinitial_x,theta0,a & 6: method to be used.
        throw runtime_error("There aren't enough input parameters.");
        return -1;
    }
    
    double initial_x=atof(argv[1]);
    double initial_y=atof(argv[2]);
    double initial_theta=atof(argv[3])*M_PI/180;
    double a=atof(argv[4]);
    string method= argv[5];
    
    if (abs(a)>1){// Unphysical black hole spin parameter.
        throw runtime_error("A black hole cannot exist with this spin value");
        return -1;
    }
    
    if (method != "Simple" && method != "Disk" && method != "NoDisk"){
        throw runtime_error("The method to be followed must be one of Simple,Disk,and NoDisk.");
        return -1;
    }
    
    // These parameters are unchanging for the Photon orbit
    // d=observer distance (approximation as in reality -> infinity)
    // F involved in time step in algorithm
    // n is a maximum number of steps per integration per photon,exists to avoid algorithm getting stuck. (Happens when the photon reaches the event horizon).
    
    double d=1000;
    double F=1.0/128.0;
    const int n=5000;
    
    // these will be updated after photon integration with the
    // points in the disk where the photon originated from
    
    double rf=0;
    double phif=0;
    double thetaf=0;
        
    example_photon.compute_orbit(a,initial_x,initial_y,initial_theta,d,F,n,rf,phif,thetaf,method); // solve the orbit equations for a single photon
    
    cout<<rf<<endl;
    cout<<thetaf<<endl;
    cout<<phif<<endl;
    return 0;
};
