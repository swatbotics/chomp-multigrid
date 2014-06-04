/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "Chomp.h"
#include "Constraint.h"
#include "ConstraintFactory.h"
#include <getopt.h>
#include <stdio.h>

#ifdef MZ_HAVE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#endif

using namespace chomp;


class CircleConstraint : public Constraint {
public:

  virtual size_t numOutputs(){
    return 1;
  }

  virtual void evaluateConstraints(const MatX& qt,
                                   MatX& h, 
                                   MatX& H){

    assert((qt.cols()==1 && qt.rows()==2) || (qt.cols() == 2 && qt.rows() == 1));

    h.conservativeResize(1,1);
    H.conservativeResize(1,2);
    
    h(0) = mydot(qt,qt) - 4; 
    
    for(int i=0; i<2; i++){
      H(i) = 2*qt(i);
    }

  }

};

class CircleFactory: public ConstraintFactory {
public:

  virtual Constraint* getConstraint(size_t t, size_t total) {

    if (4*(t+1)<(total+1) || 4*(t+1)>3*(total+1)) {
      
      return new NullConstraint();
      
    } else { 

      return new CircleConstraint();
      
    } 

  }


};

#ifdef MZ_HAVE_CAIRO

class PdfEmitter: public DebugChompObserver {
public:

  bool dumpOpt;
  int count;

  double sz, scl, m;
  cairo_surface_t* surface;
  cairo_t* cr;

  const char* filename;

  PdfEmitter(const char* f): dumpOpt(false), count(0), filename(f) {

    sz = 4*72; 
    scl = sz/10;
    m = 0;

    surface = cairo_pdf_surface_create(filename,
                                       sz+2*m, sz+2*m);

    cr = cairo_create(surface);

  }

  virtual ~PdfEmitter() {

    cairo_surface_destroy(surface);
    cairo_destroy(cr);

    std::cout << "wrote " << filename << "\n\n";
    
  }

  virtual int notify(const Chomp& chomper, 
                     ChompEventType event,
                     size_t iter,
                     double curObjective,
                     double lastObjective,
                     double hmag) {

    DebugChompObserver::notify(chomper, event, iter, 
                               curObjective, lastObjective, hmag);

    if (! ( (event == CHOMP_INIT) ||
            (event == CHOMP_FINISH) ||
            (dumpOpt) ) ) {
          
      return 0;

    }

    if (count++) { 
      cairo_show_page(cr);
    }

    double x0 = -4;
    double y0 = 6;

#define MAPX(x) (((x)-x0)*scl+m)
#define MAPY(y) ((y0-(y))*scl+m)

    cairo_set_line_width(cr, 1);
    cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
    cairo_arc(cr, MAPX(0), MAPY(0), 2*scl, 0, 2*M_PI);
    cairo_stroke(cr);

    for (int i=-1; i<=chomper.N; ++i) {
      MatX pi;
      if (i < 0) { 
        pi = chomper.q0;
      } else if (i >= chomper.N) {
        pi = chomper.q1;
      } else {
        pi = chomper.xi.row(i);
      }
      double u = double(i+1) / (chomper.N+1);
      cairo_set_source_rgb(cr, 1-u, 0, u);
      cairo_arc(cr, MAPX(pi(0)), MAPY(pi(1)), 2, 0, 2*M_PI);
      cairo_stroke(cr);
    }

    cairo_set_line_width(cr, 2);
    cairo_rectangle(cr, MAPX(-4), MAPY(6), 10*scl, 10*scl);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_stroke(cr);

      
    return 0;

  }
  

};

#endif

void generateInitialTraj(int N, MatX& xi, MatX& q0, MatX& q1) {

  xi.resize(N, 2);
  q0.resize(1, 2);
  q1.resize(1, 2);

  q0 << -3, 5;
  q1 << 5, -3;

  for (int i=0; i<N; ++i) {
    xi.row(i) = (i+1)*(q1-q0)/(N+1) + q0;
  }

}


void usage(int status) {
  std::ostream& ostr = status ? std::cerr : std::cout;
  ostr <<
    "usage: circle_demo OPTIONS\n"
    "\n"
    "OPTIONS:\n"
    "\n"
    "  -n, --num-initial        Number of steps for initial trajectory\n"
    "  -t, --num-final          Minimum number of steps for final trajectory\n"
    "  -l, --no-local           Disable local smoothing\n"
    "  -g, --no-global          Disable global smoothing\n"
    "  -m, --no-multigrid       Disable multigrid computation\n"
    "  -e, --error-tol          Relative error tolerance\n"
    "  -a, --alpha              Step size for CHOMP\n"
    "  -p, --pdf                Output PDF's\n"
    "      --help               See this message.\n";
  exit(status);
}
    

int main(int argc, char** argv) {

  int N = 15;
  int Nmax = 127;

  bool doMultigrid = true;
  bool doGlobalSmooth = true;
  bool doLocalSmooth = true;
  bool doPDF = false;
  double alpha = 0.05;
  double errorTol = 1e-7;

  const struct option long_options[] = {
    { "num-initial",       required_argument, 0, 'n' },
    { "num-final",         required_argument, 0, 't' },
    { "error-tol",         required_argument, 0, 'e' },
    { "alpha",             required_argument, 0, 'a' },
    { "no-multigrid",      no_argument,       0, 'm' },
    { "no-global",         no_argument,       0, 'g' },
    { "no-local",          no_argument,       0, 'l' },
    { "pdf",               no_argument,       0, 'p' },
    { "help",              no_argument,       0, 'h' },
    { 0,                   0,                 0,  0  }
  };

  const char* short_options = "n:t:e:a:mpglh";
  int opt, option_index;

  while ( (opt = getopt_long(argc, argv, short_options, 
                             long_options, &option_index) ) != -1 ) {

    switch (opt) {
    case 'n':
      N = atoi(optarg);
      break;
    case 't':
      Nmax = atoi(optarg);
      break;
    case 'e':
      errorTol = atof(optarg);
      break;
    case 's':
      alpha = atof(optarg);
      break;
    case 'm':
      doMultigrid = false;
      break;
    case 'g':
      doGlobalSmooth = false;
      break;
    case 'l':
      doLocalSmooth = false;
      break;
    case 'p':
      doPDF = true;
      break;
    case 'h':
      usage(0);
      break;
    default:
      usage(1);
      break;
    }

  }

  if (!doMultigrid) { 
    N = Nmax; 
    doLocalSmooth = false;
  }
  
  std::cout << "about to optimize with settings:\n";
  std::cout << "  init n:           " << N << "\n";
  std::cout << "  final n:          " << Nmax << "\n";
  std::cout << "  step size:        " << alpha << "\n";
  std::cout << "  error tol:        " << errorTol << "\n";
  std::cout << "  multigrid:        " << (doMultigrid ? "ON" : "OFF") << "\n";
  std::cout << "  local smoothing:  " << (doLocalSmooth ? "ON" : "OFF") << "\n";
  std::cout << "  global smoothing: " << (doGlobalSmooth ? "ON" : "OFF") << "\n\n";

  MatX xi, q0, q1;

  generateInitialTraj(N, xi, q0, q1);

  CircleFactory factory;

  Chomp chomper(&factory, xi, q0, q1, Nmax, alpha, errorTol);

  DebugChompObserver obs;
  chomper.observer = &obs;

#ifdef MZ_HAVE_CAIRO
  PdfEmitter* pobs = 0;

  if (doPDF) {

    char filename[1024];

    snprintf(filename, 1024, "circle_%04d_%04d_%s_%s.pdf",
             chomper.minN, chomper.maxN, 
             doLocalSmooth ? "local" : "nolocal", 
             doGlobalSmooth ? "global" : "noglobal");

    pobs = new PdfEmitter(filename);

    pobs->dumpOpt = (doMultigrid == false);
    chomper.observer = pobs;
  }
#endif

  chomper.solve(doGlobalSmooth, doLocalSmooth);

#ifdef MZ_HAVE_CAIRO
  delete pobs;
#endif

  return 0;

}
