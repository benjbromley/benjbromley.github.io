#ifndef DR2_H
#define DR2_H
#include "vec3d.h"
#include "sobseq.h"
#include "gaxgravity.h"
#include <vector>
#include <map>
#include <algorithm>


using namespace std; //fuck...

const double pi = acos(-1.0);
const double deg = pi/180.0;
template <class T> inline T sqr(const T x) { return x*x; }

const vec3d Rsungc_Marchetti(-8.2,0.0,0.025); // kpc
const vec3d Vsunlsr_Marchetti(14.0,12.24,7.25); // # km/s refs. in Marchetti
const vec3d Vlsrgc_Marchetti(0.0,238.0,0);

const vec3d Rsungc_Hattori(-8.0,0.0,0.0); // kpc
const vec3d Vsunlsr_Hattori(11.1,12.24,7.25); // # sonrich 2010
const vec3d Vlsrgc_Hattori(0,220.0,0);  // Hogg 2005; Bovy  2012; Reid 2014

const vec3d Rsungc_us(-8.0,0.0,0.0); // kpc
const vec3d Vsunlsr_us(11.1,12.24,7.25); // # sonrich 2010
const vec3d Vlsrgc_us(0,235.0,0);  // Hogg  2005; Bovy 2012;  2014

vec3d Rsungc = Rsungc_us;
vec3d Vsunlsr = Vsunlsr_us;
vec3d Vlsrgc = Vlsrgc_us;

const double kconvert=4.74047; //4.74047???

const double raG=192.85948*deg; 
const double decG=27.12825*deg;
const double lNGP=122.93192*deg;

void set_Solar_LSR(string model)
{
  if (model=="Marchetti") {
    Rsungc = Rsungc_Marchetti;
    Vsunlsr = Vsunlsr_Marchetti;
    Vlsrgc = Vlsrgc_Marchetti;
  } else if (model=="Hattori") {
    Rsungc = Rsungc_Hattori;
    Vsunlsr = Vsunlsr_Hattori;
    Vlsrgc = Vlsrgc_Hattori;
  } else if (model=="us") {
    Rsungc = Rsungc_us;
    Vsunlsr = Vsunlsr_us;
    Vlsrgc = Vlsrgc_us;
  } else if (model=="simple") {
    Rsungc = vec3d(8,0,0);
    Vsunlsr = vec3d(0,235,0);
  } else {
    cout << " bad Solar motion model " << endl;
    exit(255);
  }
  return;
}

class dr2dat
{
public:
  char designation[20];
  double ra,dec,pax,pmra,pmdec,rv;
  double ra_err,dec_err,pax_err,pmra_err,pmdec_err,rv_err;
  // double ra_dec_corr, ra_pax_corr, ra_pmra_corr, ra_pmdec_corr;
  // double dec_pax_corr, dec_pmra_corr, dec_pmdec_corr;
  double pax_pmra_corr,pax_pmdec_corr,pmra_pmdec_corr;
  double bp_rp,bp_g,magg,magbp,brxf,agof,achi2,axns,mvfa;
  int angood,vpu,rvnt;
  char dup;
};

class dr2datx : public dr2dat
{
public:
  double Gmag,Gmag_err;
  char spectype[8];
  char id[8];
 
};

inline void eq2gc(double& l, double& b, const double ra, const double dec) {
  // radians...
  b = asin(cos(dec)*cos(decG)*cos(ra-raG)+sin(dec)*sin(decG));
  l = lNGP-atan2(cos(dec)*sin(ra-raG),
		 sin(dec)*cos(decG)-cos(dec)*sin(decG)*cos(ra-raG));
  if (l>= 360.*deg) l-=360.*deg;
  if (l< 0.*deg) l+=360.*deg;
  return;
}

inline void pmeq2gc(double& pml, double& pmb, const double ra, 
		      const double dec, const double pmra, const double pmdec)
// radians here...
{
  // https://arxiv.org/pdf/1306.2945.pdf
  // radians...then whatever vel units...
  const double C1 = sin(decG)*cos(dec)-cos(decG)*sin(dec)*cos(ra-raG);
  const double C2 = cos(decG)*sin(ra-raG);
  const double icosb = 1.0/sqrt(C1*C1+C2*C2);
  pml = icosb*(C1*pmra + C2*pmdec);
  pmb = icosb*(-C2*pmra + C1*pmdec);
  return;
}

// johnson & Soderblom 1987
inline void Tpack(double(&T)[3][3]) {
  const double a=raG; const double ca=cos(a); const double sa=sin(a); 
  const double d=decG; const double cd=cos(d); const double sd=sin(d);
  const double h=lNGP; const double ch=cos(h); const double sh=sin(h);
  const double T1[3][3] = {{ca,sa,0.},{sa,-ca,0.},{0.,0.,1.}};
  const double T2[3][3] = {{-sd,0,cd},{0,-1,0},{cd,0,sd}};
  const double T3[3][3] = {{ch,sh,0},{sh,-ch,0},{0,0,1.}};
  double T12[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for (size_t i=0; i<3; i++) 
    for (size_t j=0; j<3; j++) 
      for (size_t k=0; k<3; k++) 
	T12[i][j] +=  T2[i][k]*T1[k][j];
  for (size_t i=0; i<3; i++) {
    for (size_t j=0; j<3; j++) {
      T[i][j] = 0.0;
      for (size_t k=0; k<3; k++) 
	T[i][j] += T3[i][k]*T12[k][j];
    }
  }
  return;
}


inline void eq2lsr3d(vec3d& XYZ, vec3d& UVW, const double pax, 
		     const double ra,  const double dec, const double rv, 
		     const double dmra, const double dmdec) {
  // units: mas, deg, deg
  // rectilinear coordinates, heliocentric
  static int init;
  static double T[3][3];
  if (!init++) Tpack(T);
  double ca=cos(ra*deg), sa=sin(ra*deg);
  double cd=cos(dec*deg), sd=sin(dec*deg);
  double xx=ca*cd/pax, yy=sa*cd/pax, zz=sd/pax;
  double req[3] = {xx,yy,zz};
  XYZ = vec3d(0,0,0);
  for (size_t i=0; i<3; i++) 
    for (size_t j=0; j<3; j++) 
	XYZ[i] +=  T[i][j]*req[j];
  double A[3][3] = {{ca*cd,-sa,-ca*sd},{sa*cd,ca,-sa*sd},{sd,0.0,cd}};
  double B[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for (size_t i=0; i<3; i++) 
    for (size_t j=0; j<3; j++) 
      for (size_t k=0; k<3; k++) 
	B[i][j] +=  T[i][k]*A[k][j];
  UVW = vec3d(0,0,0);
  double veq[3] = {rv,kconvert*dmra/pax,kconvert*dmdec/pax};
  for (size_t i=0; i<3; i++) 
    for (size_t j=0; j<3; j++) 
      UVW[i] += B[i][j]*veq[j];
}

void eq2gc3d(vec3d& XYZ, vec3d& UVW, const double pax, const double ra, 
	   const double dec, const double rv, const double dmra, 
	     const double dmdec) { // units km/s, mas, deg, mas/yr, wow
  eq2lsr3d(XYZ,UVW,pax,ra,dec,rv,dmra,dmdec);

  // Marchetti et al. 2018
  XYZ += Rsungc;
  UVW += Vsunlsr + Vlsrgc;  
}


inline
double vradcorrected(const double ra, const double dec, const double rv)
// expecting *degrees* for ra, dec;
{
  double l,b;
  eq2gc(l,b,ra*deg,dec*deg);
  vec3d UVWs= Vsunlsr+Vlsrgc;
  const double Us(UVWs[0]),Vs(UVWs[1]),Ws(UVWs[2]);
  const double sl(sin(l)),cl(cos(l)),sb(sin(b)),cb(cos(b));
  double vr = rv + Us*cb*cl + Vs*cb*sl + Ws*sb;
  return vr;
    
}

inline
double vtancorrected(const double ra, const double dec, const double pax, 
		     const double pmra, const double pmdec)
// expecting degrees for ra dec, mas for pax, mas/yr, returning km/s
{
  double l,b,pml,pmb;
  eq2gc(l,b,ra*deg,dec*deg);
  pmeq2gc(pml,pmb,ra*deg,dec*deg,pmra,pmdec);
  vec3d UVWs= Vsunlsr+Vlsrgc;
  const double Us(UVWs[0]),Vs(UVWs[1]),Ws(UVWs[2]);
  const double sl(sin(l)),cl(cos(l)),sb(sin(b)),cb(cos(b));
  double vt2 = sqr(kconvert*pml/pax - Us*sl + Vs*cl)+
    sqr(kconvert*pmb/pax - Us*cl*sb - Vs*sl*sb + Ws*cb);
  return sqrt(vt2);
}




class stats
{
 public:
  double mean,stddev,median,s1minus,s1plus;
  int pack_stats(vector<double>& soriginal) {
    mean = stddev = s1minus = s1plus = 0.0;
    if (soriginal.size()==0) return 0;
    vector<double> s = soriginal;
    sort(s.begin(),s.end());
    for (size_t i=0; i<s.size(); i++) mean += s[i];
    mean /= s.size();
    for (size_t i=0; i<s.size(); i++) stddev += sqr(s[i]-mean);
    stddev = sqrt(stddev/s.size());
    median=s[s.size()/2];
    s1minus = median-s[int(16./100.*s.size())];// abs errorbar siz
    s1plus = s[int(84./100.*s.size())]-median;
    s.clear();
    return 1;
  }
};

/* gaussrand dev gen. (Num. Rec.+minor change) */
inline double gaussrand()
{
  static int iset;
  static double gset;
  double fac,r,v1,v2;  
  if (iset == 0) {
    do {
      v1 = 2.*drand48()-1.;
      v2 = 2.*drand48()-1.;
      r = v1*v1+v2*v2;
    } while (r >= 1. || r == 0.);
    fac = sqrt(-2.*log(r)/r);
    gset = v1*fac;
    iset = 1;
    return(v2*fac);
  } else {
    iset = 0;
    return(gset);
  }
}

inline
void randgen(double& Y1, double& Y2, double& Y3, double& rv, const double  m1, 
	     const double m2, const double m3, const double mrv,
	     const double s11, const double s22, const double s33, 
	     const double srv2, 
	     const double s12, const double s13, const double s23, 
	     const bool need2init)
{
  static double L11=0.0, L12=0.0, L13=0.0, L22=0.0, L23=0.0, L33=0.0;
  
  if (need2init) { 
    L11=sqrt(s11); L12=s12/L11; L13=s13/L11;
    L22=sqrt(s22-L12*L12);
    L23=(s23-L13*L12)/L22;
    L33=sqrt(s33-(L13*L13+L23*L23));
  }
  double X1 = gaussrand();
  double X2 = gaussrand();
  double X3 = gaussrand();

  Y3=L13*X1+L23*X2+L33*X3; Y2=L22*X2+L12*X1; Y1=L11*X1;
  Y1+=m1; Y2+=m2; Y3+=m3;

  rv  = mrv + sqrt(srv2)*gaussrand();
  
  return;
}

inline
void gridgen(double& Y1, double& Y2, double& Y3, double& rv, const double  m1, 
	     const double m2, const double m3, const double mrv,
	     const double s11, const double s22, const double s33, 
	     const double srv2, 
	     const double s12, const double s13, const double s23, 
	     const bool need2init)
{
  static double L11=0.0, L12=0.0, L13=0.0, L22=0.0, L23=0.0, L33=0.0;
  static int is_set;
  if (need2init) { 
    L11=sqrt(s11); L12=s12/L11; L13=s13/L11;
    L22=sqrt(s22-L12*L12);
    L23=(s23-L13*L12)/L22;
    L33=sqrt(s33-(L13*L13+L23*L23));
    is_set = 0;
  }
  double XX[4];
  sobseq(XX,4,is_set);
  is_set = 1;
  double X1 = sqrt(-2*log(XX[0]))*cos(2*pi*XX[1]);
  double X2 = sqrt(-2*log(XX[0]))*sin(2*pi*XX[1]);
  double X3 = sqrt(-2*log(XX[2]))*cos(2*pi*XX[3]);
  double X4 = sqrt(-2*log(XX[2]))*sin(2*pi*XX[3]);

  Y3=L13*X1+L23*X2+L33*X3; Y2=L22*X2+L12*X1; Y1=L11*X1;
  Y1+=m1; Y2+=m2; Y3+=m3;
  rv = mrv + sqrt(srv2) * X4;
  return;
}

inline
void griddegen(double& Y1, double& Y2, double& Y3, double& rv, const double m1,
	       const double m2, const double m3, const double mrv,
	       const double s11, const double s22, const double s33, 
	       const double srv2, 
	       const double s12, const double s13, const double s23, 
	       const double rpow, const double L,   const bool need2init)
// this builds in a distance prior to get pax, so that 1/pax recovers
// the distance estimate itsel.
// so unlike prev function this one is specific to m1=pax.
{
  static double L11=0.0, L12=0.0, L13=0.0, L22=0.0, L23=0.0, L33=0.0;
  static int is_set;
  static double ddmin, ddmax;
  static vector<double> cdfd(0), dd(0);

  if (need2init) { 
    L11=sqrt(s11); L12=s12/L11; L13=s13/L11;
    L22=sqrt(s22-L12*L12);
    L23=(s23-L13*L12)/L22;
    L33=sqrt(s33-(L13*L13+L23*L23));

    ddmin = 0.001; ddmax = 250.0; int ndd = 5000;
    dd.resize(ndd); cdfd.resize(ndd);
    dd[0] = (ddmin);
    double dlogdd = log(ddmax/ddmin)/(ndd-1);
    for (int i=1; i<ndd; i++) dd[i] = exp(log(dd[0])+i*dlogdd);
    cdfd[0] = 0.0;
    for (int i=1; i<ndd; i++) {
      double r = 0.5*(dd[i-1]+dd[i]);
      double dr = (dd[i]-dd[i-1]);
      double tmp = (1/r-m1)/sqrt(s11);
      double pdfx = pow(r,rpow)*exp(-0.5*tmp*tmp - r/L);
      cdfd[i] = cdfd[i-1] + dr * pdfx;
    }
    double normfac = cdfd.back();
    for (int i=0; i<ndd; i++) cdfd[i] /= normfac;
    is_set = 0;
  }
  double XX[4];
  sobseq(XX,4,is_set);
  is_set = 1;
  double X1 = sqrt(-2*log(XX[0]))*cos(2*pi*XX[1]);
  double X2 = sqrt(-2*log(XX[0]))*sin(2*pi*XX[1]);
  double X3 = sqrt(-2*log(XX[2]))*cos(2*pi*XX[3]);
  double X4 = sqrt(-2*log(XX[2]))*sin(2*pi*XX[3]);
  double U1 = 0.5*(erf(X1/sqrt(2.))+1.0);
  // invert d here:
  vector<double>::iterator it = upper_bound(cdfd.begin(),cdfd.end(),U1);
  if (it==cdfd.end()) --it;
  size_t ii = distance(cdfd.begin(),it);
  double w = (U1-cdfd[ii-1])/(cdfd[ii]-cdfd[ii-1]);
  double d = (1.0-w)*dd[ii-1]+w*dd[ii];
  X1 = (1.0/d-m1)/L11;
  Y3=L13*X1+L23*X2+L33*X3; Y2=L22*X2+L12*X1; Y1=L11*X1;
  Y1+=m1; Y2+=m2; Y3+=m3;
  rv = mrv + sqrt(srv2) * X4;
  return;
}


inline double distanceinfer(const double paxtry, const double paxbar, 
			    const double paxerr, const double rpow, 
			    const double L, const bool do_initialization)
{
  static double ddmin, ddmax, ucut;
  static vector<double> cdfd(0), dd(0);
  if (do_initialization) {
    ddmin = 0.001; ddmax = 250.0; int ndd = 5000;
    dd.resize(ndd);
    cdfd.resize(ndd);
    dd[0] = (ddmin);
    double dlogdd = log(ddmax/ddmin)/(ndd-1);
    for (int i=1; i<ndd; i++) dd[i] = exp(log(dd[0])+i*dlogdd);
    //dd[0] + exp(i*dlogdd);
    // for (int i=1; i<ndd; i++) dd[i] = dd[0] + i*(ddmax-ddmin)/(ndd+1);
    cdfd[0] = 0.0;
    for (int i=1; i<ndd; i++) {
      double r = 0.5*(dd[i-1]+dd[i]);
      double dr = (dd[i]-dd[i-1]);
      double tmp = (1/r-paxbar)/paxerr;
      double pdfx = pow(r,rpow)*exp(-0.5*tmp*tmp - r/L);
      cdfd[i] = cdfd[i-1] + dr * pdfx;
    }
    double normfac = cdfd.back();
    for (int i=0; i<ndd; i++) cdfd[i] /= normfac;
    ucut = 0.5*(erf((-paxbar)/sqrt(2.0)/paxerr)+1.0);  
  }
  if (paxtry<=0) return -1e99;
  double u = 0.5*(erf((paxtry-paxbar)/sqrt(2.0)/paxerr)+1.0);
  //  cout << " u raw " << u << " " << paxtry << " " << paxerr << endl;
  //cout << " xxx raw " << (paxtry-paxbar)/sqrt(2.0)/paxerr << endl;
  //cout << erf(-1.0) << " " << erf(0.0);
  u -= ucut; u /= (1.0-ucut); u = 1-u;
  vector<double>::iterator it = upper_bound(cdfd.begin(),cdfd.end(),u);
  if (it==cdfd.end()) return -1e99;
  size_t ii = distance(cdfd.begin(),it);
  double w;
  w = (u-cdfd[ii-1])/(cdfd[ii]-cdfd[ii-1]);
  double d = (1.0-w)*dd[ii-1]+w*dd[ii];
  return d;
}


inline int getvesc(vector<double>&d, vector<double>&v, string fin)
{
  d.clear(); v.clear();
  char line[BUFSIZ];
  FILE *fp = fopen(fin.c_str(),"r");
  if (fp == NULL) { cout << " unable to open "<< fin << endl; exit (255); }
  while (fgets(line,BUFSIZ,fp)) {
    if (line[0]=='#') continue;
    double x,y;
    int j = sscanf(line,"%lf %lf",&x,&y);
    if (j != 2) { cout << " read err " << fin << endl; exit(255); }
    x /= 1e3; // pc to kpc
    y /= 1e5; // cm/s to km/s
    d.push_back(x);
    v.push_back(y);
  }
  return (int)d.size();
}


double vesc_interp(const double d, vector<double> dd, vector<double> vv,
		   const int need2init)
{
  static map<double, size_t> dtoi;
  if (need2init) { 
    for (size_t j=0; j<dd.size(); j++) {
      dtoi[dd[j]] = j; 
    }
  }
  if (d<=dd[0]) return vv.front();
  if (d>=dd.back()) return 0.0;
  size_t i = dtoi.lower_bound(d)->second;
  double w = (dd[i]-d)/(dd[i]-dd[i-1]);
  // cout << i << " " << d << " " << dd[i-1] << " " << dd[i] <<  endl;
  if (fabs(w)>1.001)  { cout << " prob in interp: w=" << w << endl; exit(255);}
  return w*vv[i-1]+(1.0-w)*vv[i];
}

static string _fesc = "esc.dat";
static vector<double> _Desc(0),_Vesc(0);

double vesc_gaxbdh_k(const vec3d R) // kpc, km/s
{
  using namespace gaxgravity;
  if (R.norm()>=250.0) return 0.0;
  if (R.norm()==0.0) return vesc_gaxbdh_k(vec3d(0,0,0.0001));
  double vesc = sqrt(2.0*(Phigax_k(vec3d(0,0,250*kpc))-Phigax_k(R*kpc)));
  return vesc/1e5;
}

double vesc_from_fesc(const vec3d R) { // R in kpc, returns in km/s
  if (_Desc.size()==0) {
    (void)getvesc(_Desc,_Vesc,_fesc);
    (void)vesc_interp(8.0,_Desc,_Vesc,1); // initialize map in look-up table.
  }
  return vesc_interp(R.norm(),_Desc,_Vesc,0);
}

inline void teststats()
{
  vector<double> x(0),y(0),z(0);
  const int n=1000;
  for (int i=0; i<n; i++) {
    x.push_back(gaussrand());
  }
  stats xs,ys,zs; 
  xs.pack_stats(x);
  cout << xs.mean << " " << xs.stddev << " " << xs.median << " " <<
    xs.s1minus << " " << xs.s1plus << endl;
  x.clear();
  for (int i=0; i<n; i++) {
    double a,b,c,d;
    randgen(a,b,c,d,1.0,2.0,3.0,4.0,1.0,4.0,9.0,16.0,0.18,0.13,0.23,!i);
    //    gridgen(a,b,c,1.0,2.0,3.0,1.0,4.0,9.0,0.18,0.13,0.23,!i);
    x.push_back(a);
    y.push_back(b);
    z.push_back(c);
  }
  xs.pack_stats(x);
  ys.pack_stats(y);
  zs.pack_stats(z);
  double sxy = 0;
  for (int i=0; i<n; i++) {
    sxy += (x[i]-xs.mean)*(y[i]-ys.mean);
  }
  sxy/=n;
  cout << xs.mean << " " << xs.stddev << " " << xs.median << " " <<
    xs.s1minus << " " << xs.s1plus << endl;
  cout << ys.mean << " " << ys.stddev << " " << ys.median << " " <<
    ys.s1minus << " " << ys.s1plus << endl;
  cout << zs.mean << " " << zs.stddev << " " << zs.median << " " <<
    zs.s1minus << " " << zs.s1plus << endl;
  cout << sxy << endl; 
  exit(0);
}

const int Ngcsamp = 4*1024;
bool use_sobseq = true;
bool use_prior = false;
double hyperparm_L = 1e99;
double hyperparm_rpow = -2.0; // this is no parm...

inline void gcstats(stats& ds,stats& rs, stats& vs, stats& as, stats& vrs,
		    stats& vts, double& prob_ub, double& prob_ubvr, 
		    double& prob_ubvt, dr2dat& p, 
		    double (*vescfunc)(const vec3d R), const int N)
{
  //  const int N = 4*1024;
  vector<double> dh(0),r(0),v(0),a(0),vr(0),vt(0);
  int nunbound(0),nunboundvr(0),nunboundvt(0);
  // bool initprior = true;
  for (int i=0; i<N; i++) { 
    double pax,pmra,pmdec,rv;
    if (use_sobseq==false) {
      randgen(pax,pmra,pmdec,rv,p.pax,p.pmra,p.pmdec,p.rv,sqr(p.pax_err),
	      sqr(p.pmra_err), sqr(p.pmdec_err), sqr(p.rv_err),
	      p.pax_pmra_corr*p.pax_err*p.pmra_err,
	      p.pax_pmdec_corr*p.pax_err*p.pmdec_err,
	      p.pmra_pmdec_corr*p.pmra_err*p.pmdec_err,!i);
    } else {
      if (use_prior) {
	griddegen(pax,pmra,pmdec,rv,p.pax,p.pmra,p.pmdec,p.rv,sqr(p.pax_err),
		sqr(p.pmra_err),sqr(p.pmdec_err), sqr(p.rv_err), // vars...
		p.pax_pmra_corr*p.pax_err*p.pmra_err,
		p.pax_pmdec_corr*p.pax_err*p.pmdec_err,
		p.pmra_pmdec_corr*p.pmra_err*p.pmdec_err,
		  hyperparm_rpow,hyperparm_L,!i);
      } else {
	gridgen(pax,pmra,pmdec,rv,p.pax,p.pmra,p.pmdec,p.rv,sqr(p.pax_err),
		sqr(p.pmra_err),sqr(p.pmdec_err), sqr(p.rv_err), // vars...
		p.pax_pmra_corr*p.pax_err*p.pmra_err,
		p.pax_pmdec_corr*p.pax_err*p.pmdec_err,
		p.pmra_pmdec_corr*p.pmra_err*p.pmdec_err,!i);
      }
    }
    vec3d R,V;
    if (pax>0.0) {
      // old: now griddegen returns a pax that is 1/d, where d 
      // is drawn from posterior...
      // if (use_prior) {
      //   double d = distanceinfer(pax,p.pax,p.pax_err,hyperparm_rpow,
      //			 hyperparm_L,initprior);
      //   initprior = false;
      //   pax = 1.0/d;
      // }
      eq2gc3d(R,V,pax,p.ra,p.dec,rv,pmra,pmdec);
      dh.push_back(1./pax); // will have been adjusted....
      r.push_back(R.norm());
      v.push_back(V.norm());
      a.push_back(acos(R.unit()*V.unit()));
      double Vr = vradcorrected(p.ra,p.dec,rv);
      double Vt = vtancorrected(p.ra,p.dec,pax,pmra,pmdec);
      vr.push_back(Vr);
      vt.push_back(Vt);
      double vescR = vescfunc(R);
      if (V.norm() > vescR) nunbound++;
      if (fabs(Vr) > vescR) nunboundvr++;
      if (Vt > vescR) nunboundvt++;
    }
  }
  ds.pack_stats(dh);
  rs.pack_stats(r);
  vs.pack_stats(v);
  as.pack_stats(a);
  vrs.pack_stats(vr);
  vts.pack_stats(vt);
  prob_ub = float(nunbound)/r.size();
  prob_ubvr = float(nunboundvr)/r.size();
  prob_ubvt = float(nunboundvt)/r.size();
  v.clear();
  r.clear();
  return;
}

// new stuff for streamlined searches...
// should put this in dr2reads...



inline
int astrometry_flag(dr2dat& p)
// check_phometry: 0 => don't bother; 1=>gmag only, 2=> gmag and bp-rp
// note: gmag and color can influence the quality of the astronometry
{
  int ret = 0;
  if (isnan(p.pax) && isnan(p.agof) && isnan(p.axns) && isnan(p.mvfa) &&
      isnan(p.magg) && isnan(p.bp_rp) && isnan(p.achi2) && isnan(p.brxf))
    return (1<<10)-1;

  if (isnan(p.agof)) ret += 1;
  else if (p.agof >= 3.0) ret += 1;  // marchetti

  if (isnan(p.axns)) ret += 2;
  else if (p.axns > 2.0) ret += 2;

  if (isnan(p.axns) || isnan(p.mvfa)) ret += 4;
  else if (p.mvfa <= -0.23 or p.mvfa > 0.32) ret += 4;
  // if (p.vpu <= 8) return false;  // marchetti

  const int VPUmin = 10;
  if (isnan(p.vpu)) ret += 8;
  else if (p.vpu < VPUmin) ret += 8;
  // lindegren et al. 2018

  if (isnan(p.magg) || isnan(p.bp_rp) || isnan(p.achi2) || isnan(p.brxf)
      || isnan(p.angood)) {
    ret += 16 + 32;
  }  else {  
    double u = (p.achi2/(p.angood-5));
    if (u >= 1.2*max(1.0,exp(-0.2*(p.magg-19.5)))) ret += 16;
    const double br2=p.bp_rp*p.bp_rp;
    if (p.brxf<=1.0+0.015*br2 || p.brxf>=1.3+0.06*br2) ret += 32;
  }
  return ret;
}


class dr2redux {
public:
  double ra,dec,pax,paxerr;
  double pmra,pmraerr,pmdec,pmdecerr,pmrapmdeccorr;
  double bprp,gmag;
  unsigned long long des;
  int flag;
  dr2redux() {}
  
  dr2redux(dr2dat x) {
    ra = x.ra; dec = x.dec;
    pax = x.pax; paxerr = x.pax_err;
    pmra = x.pmra; pmraerr = x.pmra_err;
    pmdec = x.pmdec; pmdecerr = x.pmdec_err;
    pmrapmdeccorr = x.pmra_pmdec_corr;
    bprp = x.bp_rp; gmag = x.magg;
    char **tmp = NULL;
    des = strtoull(x.designation,tmp,10);
    flag = astrometry_flag(x);
  }
};


#endif
