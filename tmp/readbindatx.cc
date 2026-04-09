#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <dirent.h>
#include <list>
#include "dr2.h"

using namespace std;

int main(int argc, char *argv[])
{
  int select_sky = 0;
  double ra = 0.0;
  double dec = 0.0;
  double rad = 0.0;
  double ralo = 0.0;
  double rahi = 360.0;
  double declo = -90.0;
  double dechi = 90.0;
  double blo = 0.0;
  double bhi = 360.0;
  double llo = -90.0;
  double lhi = 90.0;
  
  bool verbose = false;
  string fin(""),fout("");

  if (argc == 1) {
    cout <<  "usage: cmd infile options..." << endl;
    exit(255);
  } else {
    while (argc > 1) {
      string str = *++argv;
      if (str == "-verbose") {
	verbose = 1;
	--argc;
	argc -= 2;
      } else if (str == "-circ" && argc > 4) {  // -circ RA DEC RADIUS (degrees)
	select_sky += 1;
	ra = atof(*++argv); // degrees
	dec = atof(*++argv);
	rad = atof(*++argv);
	argc -= 4;
      } else if (str == "-boxradec" && argc>5) { // -boxradec RA_LO DEC_LO RA_HI DEC_HI
	select_sky += 2;
	ralo = atof(*++argv); // degrees
	declo = atof(*++argv);
	rahi = atof(*++argv); // degrees
	dechi = atof(*++argv);
	argc -= 5;
      } else if (str == "-boxgax" && argc > 5) { // -boxgax L_LO B_LO L_HI B_HI
	select_sky += 4;
	llo = atof(*++argv); // degrees
	blo = atof(*++argv);
	lhi = atof(*++argv); // degrees
	bhi = atof(*++argv);
	argc -= 5;
      } else if (str == "-fout" && argc > 2) {
	fout = (*++argv);
	argc -= 2;
      } else if (fin == "" && str[0] != '-') {
	fin = str;
	argc--;
      } else {
	cout << " bad command line argument? [" << str << "]" << endl;
	exit(255);
      }
    }
  }

  int nread = 0;
  int nwrite = 0;

  bool binout = false;
  FILE *fo = NULL;
  if (fout.size()>0) {
    fo = fopen(fout.c_str(),"w");
    binout = true;
  }
  
  FILE *fi = fopen(fin.c_str(),"r");
  if (fi == NULL) { 
    cout << " file read failure..." << fin << endl;
    exit(255);
  }
  cout << "# opened file " << fin << endl;

  //  vector<dr2redux> b(0);
  {
    int printhdr = 0;
    double dr(dec*deg),ar(ra*deg);
    double xr(cos(ar)*cos(dr)),yr(sin(ar)*cos(dr)),zr(sin(dr));
    vec3d ez = vec3d(xr,yr,zr).unit();
    vec3d ex = vec3d(-yr,xr,0).unit();
    vec3d ey = ez % ex;

    dr2redux bx;
    while (fread(&bx,sizeof(dr2redux),1,fi)==1) {
      nread++;
      if (select_sky & 1) {
	double d1(bx.dec*deg),d2(dec*deg);
	double a1(bx.ra*deg),a2(ra*deg);
	double cosdel=sin(d1)*sin(d2)+cos(d1)*cos(d2)*cos(a1-a2);
	if (cosdel<cos(rad*deg)) continue;
      }
      if (select_sky & 2) {
	if (!(bx.dec>=declo && bx.dec<=dechi && bx.ra>=ralo && bx.ra<=rahi)) continue;
      }
      if (select_sky & 4) {
	double ll,bb; eq2gc(ll,bb,bx.ra*deg,bx.dec*deg);
	ll /= deg; bb /= deg;
	if (llo>lhi) { 
	  if (!(ll>=lhi && ll<=llo && bb>=blo && bb<=bhi)) continue;
	} else {
	  if (!(ll>=llo && ll<=lhi && bb>=blo && bb<=bhi)) continue;
	}
      }
      if (false) { // junk code....
	double d(bx.dec*deg),a(bx.ra*deg);
	double x(cos(a)*cos(d)),y(sin(a)*cos(d)),z(sin(d));
	vec3d del(x-xr,y-yr,z-zr);
      }
      // double dx(del*ex),dy(del*ey);

      if (binout) {
	fwrite((void*)&bx,sizeof(bx),1,fo);
      } else {
	if (!printhdr) {
	  printf("# ra dec pax paxerr pmra pmraerr pmdec pmdecerr ");
	  printf("pmrapmdeccorr bprp gmag flag des\n");
	}
	if (printhdr == 0) { 
	  printhdr = 1;
	}
	printf("%17.15g %17.15g ",bx.ra,  bx.dec);
	printf("%10.8g %10.8g ",bx.pax,   bx.paxerr);
	printf("%10.8g %10.8g ",bx.pmra,  bx.pmraerr);
	printf("%10.8g %10.8g ",bx.pmdec, bx.pmdecerr);
	printf("%10.8g  ",bx.pmrapmdeccorr);
	printf("%10.8g %10.8g %4d ",bx.bprp,bx.gmag,bx.flag);
	// printf("%10.8g %10.8g ",dx,dy);
	printf("%llu ",bx.des);
	printf("\n");
      // fflush(stdout);
      }
      nwrite++;
      // b.push_back(bx);
    }
  }
  fclose(fi);
  if (binout) fclose(fo);
  cout << "# read: " << nread << endl;
  cout << "# kept: " << nwrite << endl;
  return 0;
}

