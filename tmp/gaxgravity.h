#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vec3d.h"


namespace gaxgravity {
  
  using namespace std;

  const double MSolar = 1.98855e+33;
  const double GNewt = 6.67408e-8;
  const double pc = 3.0857e18;
  const double kpc = pc*1e3;
  const double pi = acos(-1.0);

  template <class T> inline T sqr(const T x) { return x*x; }

  vec3d accgax_b(vec3d r, double& phi);
  vec3d accgax_k(vec3d r, double& phi);

  // galaxy mass model params, bkg=bromley,kenyon,geller,et al.
  const double Cbkg = 13962.0*MSolar/(pc*pc*pc);
  const double rcbkg = 8*pc;

  vec3d accgax_b(vec3d r, double& phi)
  {
    // bkg model:
    const double rc2bkg = rcbkg*rcbkg;
    double rr = r.norm();
    phi=2*pi*GNewt*Cbkg*rc2bkg*(2*rcbkg*atan(rr/rcbkg)/rr+log(1+rr*rr/rc2bkg));
    return -(4*pi*GNewt*Cbkg*rc2bkg*(rr-rcbkg*atan(rr/rcbkg))/(r*r))*r.unit();
  }

  // the magnificent seven from KBGB 14:

  const double Mbulge=3.75e9*MSolar;
  const double abulge = 0.105*kpc;
  const double Mdisk = 6e10*MSolar; //
  const double adisk = 2.75*kpc;
  const double bdisk = 0.3*kpc;
  const double Mhalo = 1e12*MSolar;
  const double ahalo = 20*kpc;
  const double G = GNewt;
  const double Mbh = 3.5e6*MSolar;

  double Phib(const double r)
  {
    return -G*Mbulge/(r+abulge);
  }

  vec3d accb(vec3d& r)
  {
    double rn = r.norm();
    return -G*Mbulge*r/(rn*sqr(rn+abulge));
  }

  double Phid(const double R, const double z)
  {
    return -G*Mdisk/sqrt(R*R+sqr(adisk+sqrt(z*z+bdisk*bdisk)));
  }

  vec3d accd(vec3d& r)
  {
    vec3d rp(r.x,r.y,0);
    vec3d z(0,0,r.z);
    double R2 = r.x*r.x+r.y*r.y;
    double z2 = r.z*r.z;
    double rootz2b2 = sqrt(z2+bdisk*bdisk);
    double Razd = (R2+sqr(adisk+rootz2b2));
    double Rdp = (Razd*sqrt(Razd));
    double ax = r.x/Rdp;
    double ay = r.y/Rdp;
    double az = r.z*(adisk+rootz2b2)/Rdp/rootz2b2;
    return -G*Mdisk*vec3d(ax,ay,az);
  }


  double Phih(const double r)
  {
    return -G*Mhalo*log(1+r/ahalo)/r;
  }

  vec3d acch(vec3d& r)
  {
    const double rs = ahalo;
    double rn = r.norm();
    return -G*Mhalo*r/(rn*rn*rn)*(log(1+rn/rs)-rn/(rs+rn));
  }

  vec3d accgax_k(vec3d r, double& phi)
  {
    double rn = r.norm();
    double R = sqrt(r.x*r.x+r.y*r.y);
    phi = Phih(rn)+Phid(R,r.z)+Phib(rn);
    phi += -G*Mbh/rn;
    return acch(r)+accb(r)+accd(r) - G*Mbh*r/(rn*rn*rn);
  }

  double Phigax_k(vec3d r)
  {
    double rn = r.norm();
    double R = sqrt(r.x*r.x+r.y*r.y);
    double phi = Phih(rn)+Phid(R,r.z)+Phib(rn);
    phi += -G*Mbh/rn;
    return phi;
  }
}
