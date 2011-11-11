#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

/*
 * CharmBH: Parameters.h
 * File containing structure to hold simulation parameters,
 * and the means to extract these parameters from command
 * line arguments.
 */

#include "defines.h"
#include "charm++.h"

#include <iostream>
#include <stdio.h>
#include <string>
#include <map>

#include "pup_stl.h"

using namespace std;

struct Parameters {
  // Input file name
  char *filename;
  int nchars;

  // Gravitational force calculation parameters
  Real theta;
  Real dtime;
  Real dthf;
  Real epssq;
  Real tolsq;

  // Parallel simulation parameters (see Main.cpp for explanations)
  int numTreePieces;
  int numParticles;
  int ppc;
  int ppb;

  int yieldPeriod;
  int cacheLineSize;

  int iterations;
  int balancePeriod;

  int combineFlushCount;
  int combineFlushPeriod;

  bool doPrintAccel;

  //int branchFactor;

  /*
   * Define a pup() method so that this structure can be 
   * declared a read-only, and transferred to all PEs from
   * the PE of the mainchare after initialization.
   */
  void pup(PUP::er &p){
    p | nchars;
    if(p.isUnpacking()){
      filename = new char[nchars+1];
    }
    PUParray(p,filename,nchars+1);

    p | numTreePieces;
    p | numParticles;
    p | dtime;
    p | dthf;
    p | tolsq;
    p | epssq;
    //p | branchFactor;
    p | ppc;
    p | ppb;
    p | yieldPeriod;
    p | theta;
    p | cacheLineSize;
    p | iterations;
    p | balancePeriod;
    p | doPrintAccel;

    p | combineFlushCount;
    p | combineFlushPeriod;
  }

  /*
   * Methods to extract parameters from command line arguments. There is 
   * one method for each primitive type of argument (int, bool, Real, etc.)
   * These have been adapted from the SPLASH "barnes" benchmark.
   */
  void extractParameters(int argc, char **argv, map<string,string> &tab){
    for(int i = 1; i < argc; i++){
      string arg = string(argv[i]);
      size_t mpos = arg.find("-");
      size_t pos = arg.find("=");
      size_t len = arg.length();
      if(mpos != 0) continue;

      string key, val; 
      if(pos != string::npos){
        key = arg.substr(mpos+1,pos-mpos-1);
        val = arg.substr(pos+1,len-pos-1);
      }
      else{
        key = arg.substr(mpos+1,len-mpos-1); 
        val = "";
      }
      CkPrintf("arg %s mpos %d pos %d len %d key %s val %s\n", arg.c_str(), mpos, pos, len, key.c_str(), val.c_str());

      tab[key] = val;
    }
  }

  string getparam(string name, map<string,string> &table)
  {
    map<string,string>::iterator it = table.find(name);
    if(it != table.end()){
      return it->second;
    }
    //cerr << "getparam: " << name.c_str() << " unknown" << endl;
    return string();
  }

  /*
   * GETIPARAM, ..., GETDPARAM: get int, long, bool, or double parameters.
   */

  int getiparam(string name, int def, map<string,string> &tab)
  {
    string val;

    val = getparam(name,tab);
    if(val.empty())
      return def;
    else
      return (atoi(val.c_str()));
  }

  long getlparam(string name, map<string,string> &tab)
  {
    string val;

    val = getparam(name,tab);
    if(val.empty())
      return -1;
    else 
      return (atol(val.c_str()));
  }

  bool getbparam(string name, map<string,string> &tab)
  {
    string val;

    val = getparam(name,tab);
    if (strchr("tTyY1", *(val.c_str())) != 0) {
      return (true);
    }
    if (strchr("fFnN0", *(val.c_str())) != 0) {
      return (false);
    }
    fprintf(stderr,"getbparam: %s=%s not bool\n", name.c_str(), val.c_str());
    return false;
  }

  Real getrparam(string name, Real default_value, map<string,string> &tab)
  {
    string val;

    val = getparam(name,tab);
    if(val.empty())
      return default_value;
    else 
      return (atof(val.c_str()));
  }

  /*
   * EXTRVALUE: extract value from name=value string.
   */

  string getsparam(string arg, map<string,string> &tab)
  {
    return getparam(arg,tab);
  }


};

#endif
