module tree;

import std.algorithm.iteration;
import std.stdio;
import std.range;
import std.conv;
import std.math;

import parameters;
import array2d;
import state;

immutable double pii = 3.1415926535;

class Tree
{

 public:

  int nr;
  int ny;

  double fluxin_c;
  double fluxside_c;
  double fluxtop_a_c;
  double fluxtop_d_c;

  double[] crossect;
  double[] volume;
  double[] volumeair;
  double[] volumewat;
  double[] ker1;
  double[] ker2;

  Matrix axdiff;
  Matrix eqwat;
  Matrix towater;

  State s;

  @property double time() { return s.time; }

  this(in ref Parameters p){

    // integers
    this.nr = p.nr;
    this.ny = p.ny;

    // doubles
    this.fluxin_c = 0.0;
    this.fluxside_c = 0.0;
    this.fluxtop_a_c = 0.0;
    this.fluxtop_d_c = 0.0;

    // vectors
    this.crossect.length = p.nr;
    this.volume.length = p.nr;
    this.volumeair.length = p.nr;
    this.volumewat.length = p.nr;
    this.ker1.length = p.nr;
    this.ker2.length = p.nr;

    // matrices
    this.axdiff.initialise(p.ny,p.nr);
    this.eqwat.initialise(p.ny,p.nr);
    this.towater.initialise(p.ny,p.nr);

    // state
    this.s = State(p.ny,p.nr);

    for(auto j = 0; j < p.nr; j++){
      crossect[j] = crossect_area(j,p.dr);
    }
    
    volume[0]    = p.dy * crossect[0];
    volumeair[0] = p.fractair * volume[0];
    volumewat[0] = p.fractwat * volume[0];
    ker1[0]      = double.max;
    ker2[0]      = log(2.0*p.dr/p.dr);

    // check indexing !! (was 1 based)
    double drsqr;
    for(auto j = 1; j < p.nr; j++){
      //drsqr = (to!double(j+1) * p.dr)*(to!double(j+1) * p.dr);
      //debug(3) writeln("j ",j," drsqr ",drsqr);
      //crossect[j] = (pii * drsqr) - crossect[j-1];
      volume[j]    = p.dy * crossect[j];
      volumeair[j] = p.fractair * volume[j];
      volumewat[j] = p.fractwat * volume[j];
      ker1[j] = log( (p.dr * to!double(j + 1)) /
		     (p.dr * to!double(j)) );
      ker2[j] = log( (p.dr * to!double(j + 2)) /
		     (p.dr * to!double(j + 1)) );
    }

    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
        s.cair[i,j] = p.camb;
        s.cwat[i,j] = s.cair[i,j]*p.henry;
        s.nair[i,j] = s.cair[i,j]*volumeair[j];
        s.nwat[i,j] = s.cwat[i,j]*volumewat[j];
        axdiff[i,j] = 0.0;
      }
    }

    debug(2){
      writeln("initial tree state:");
      write("crossect:");
      for(auto i = 0; i < p.nr; i++){ writef(" %s", crossect[i]); }
      writeln("");
      writeln("crossect sum: ",sum(crossect));
      write("volume:");
      for(auto i = 0; i < p.nr; i++){ writef(" %s", volume[i]); }
      writeln("");
      writeln("volume sum: ",sum(volume));
      write("ker1:");
      for(auto i = 0; i < p.nr; i++){ writef(" %s", ker1[i]); }
      writeln("");
      write("ker2:");
      for(auto i = 0; i < p.nr; i++){ writef(" %s", ker2[i]); }
      writeln("");
    }
  }

  /**
   *  used for setting areas in the constructor
   */
  pure double crossect_area(in int idx, in double dr){
    if(idx == 0){
      return pii*dr*dr;
    }else{
      return pii*pow(to!double(idx+1)*dr,2) - pii*pow(to!double(idx)*dr,2);
    }
  }

  void euler_step(in ref Parameters p){
    
    // take res steps at a time
    for(auto step = 0; step < p.res; step++){
    
      // part 1 -- radial diffusion
      for(auto i = 0; i < p.ny; i++){
	for(auto j = 0; j < p.nr; j++){

	  // set qair, radial diffusion
	  if(j == 0){
	    // innermost radial element
	    s.qair[i,j] = -2.0*pii*p.diff_r*(s.cair[i,j] - s.cair[i,j+1]) / ker2[j];
	  }else if(j == (nr - 1)){
	    // outermost radial element
	    //if(i < p.nys){
	    // root compartment, no outwards diffusion
	    //s.qair[i,j] = 2.0*pii*diff_r*(s.cair[i,j-1] - s.cair[i,j]) / ker1[j];
	    //}else{
	    // stem compartment
	    s.qair[i,j] = (2.0*pii*p.diff_r*(s.cair[i,j-1] - s.cair[i,j]) / ker1[j]) -
	      (2.0*pii*p.diff_b*(s.cair[i,j] - p.camb) / ker2[j]);
	    //}
	  }else{
	    s.qair[i,j] = (2.0*pii*p.diff_r*(s.cair[i,j-1] - s.cair[i,j]) / ker1[j]) -
	      (2.0*pii*p.diff_r*(s.cair[i,j] - s.cair[i,j+1]) / ker2[j]);
	  }

	  // set fluxdiff
	  s.fluxside[i] = 0.0;
	  //if(i > p.nys){
	  if(j == (nr - 1)){
	    s.fluxside[i] = 2.0*pii*p.diff_b*(s.cair[i,j] - p.camb) / ker2[j];
	  }
	  //}

	  // set qadv, axial advection
	  // modified: multiply every term with crossect
	  if(j < p.nrs){
	    // zero inside the duramen
	    s.qadv[i,j] = 0.0;
	  }else{
	    if(i == 0){
	      // bottom element layer receives CH4 in water
	      // s.qadv[i,j] = (p.vel*p.csoil) - (p.vel*s.cwat[i,j]);
	      // new version with element area:
	      s.qadv[i,j] = (p.vel*p.fractwat*p.csoil*crossect[j]) -
		(p.vel*s.cwat[i,j]*crossect[j]);
	      // todo: explain what happens with just this:
	      // s.qadv[i,j] = -p.vel*s.cwat[i,j];
	    }else if(i == (ny - 1)){
	      // the topmost elements simply lose CH4 in water
	      // todo: check this!
	      s.qadv[i,j] = p.vel*s.cwat[i-1,j]*crossect[j] - p.vel*s.cwat[i,j]*crossect[j];
	    }else{
	      s.qadv[i,j] = p.vel*s.cwat[i-1,j]*crossect[j] - p.vel*s.cwat[i,j]*crossect[j];
	    }
	  }	  
	  
	  // set fluxadv
	  s.fluxtop_a[j] = p.vel * s.cwat[ny-1,j] * crossect[j];

	}
      }

      // part 2 -- aksiaalinen diffusio
      for(auto i = 0; i < p.ny; i++){
	for(auto j = 0; j < p.nr; j++){
	  if(i == 0){
	    // without axial diffusion in
	    axdiff[i,j] = -(p.diff_a*crossect[j]*(s.cair[i,j] - s.cair[i+1,j]) / p.dy);
	    // with axial diffusion in
	    /*
	    axdiff[i,j] = (p.diff_a*crossect[j]*(p.csoil - s.cair[i,j]) / p.dy) -
	      (p.diff_a*crossect[j]*(s.cair[i,j] - s.cair[i+1,j]) / p.dy);
	    */
	  }else if(i == (ny - 1)){
	    axdiff[i,j] = (p.diff_a*crossect[j]*(s.cair[i-1,j] - s.cair[i,j]) / p.dy) -
	      (p.diff_a*crossect[j]*(s.cair[i,j] - p.camb) / p.dy);
	    s.fluxtop_d[j] = (p.diff_a*crossect[j]*(s.cair[i,j] - p.camb) / p.dy);
	  }else{
	    axdiff[i,j] = (p.diff_a*crossect[j]*(s.cair[i-1,j] - s.cair[i,j]) / p.dy) -
	      (p.diff_a*crossect[j]*(s.cair[i,j] - s.cair[i+1,j]) / p.dy);
	  }
	}
      }

      // part k -- kumulantit
      // do this only if flag true
      fluxin_c += p.dt * fluxin_sum(p);
      fluxside_c += p.dt * fluxside_sum();
      fluxtop_d_c += p.dt * fluxtop_d_sum();
      fluxtop_a_c += p.dt * fluxtop_a_sum();

      // part 3 -- euler first order explicit method
      for(auto i = 0; i < p.ny; i++){
	for(auto j = 0; j < p.nr; j++){
	  s.nair[i,j] = s.nair[i,j] + p.dt*(s.qair[i,j] + axdiff[i,j]);
	  s.nwat[i,j] = s.nwat[i,j] + p.dt*s.qadv[i,j];
	  s.cair[i,j] = s.nair[i,j]/volumeair[j];
	  s.cwat[i,j] = s.nwat[i,j]/volumewat[j];
	}
      }

      // part 4 -- equilibrium with dissolved CH4
      for(auto i = 0; i < p.ny; i++){
	for(auto j = 0; j < p.nr; j++){
	  eqwat[i,j] = p.henry*s.cair[i,j];
	  towater[i,j] = (eqwat[i,j] - s.cwat[i,j])*volumewat[j]*p.gamma;
	  s.nair[i,j] = s.nair[i,j] - p.dt*towater[i,j];
	  s.nwat[i,j] = s.nwat[i,j] + p.dt*towater[i,j];
	  s.cair[i,j] = s.nair[i,j]/volumeair[j];
	  s.cwat[i,j] = s.nwat[i,j]/volumewat[j];
	}
      }

      // increment time in the state object
      s.time += p.dt;

    }
  }

  pure double storage_sum(){
    double sum = 0.0;
    for(auto i = 0; i < ny; i++){
      for(auto j = 0; j < nr; j++){
	sum += s.nwat[i,j] + s.nair[i,j];
      }
    }
    return sum;
  }

  pure double[] fluxin(in ref Parameters p){
    double[] fin;
    fin.length = nr;
    for(auto j = 0; j < nr; j++){
      fin[j] = p.vel*p.fractwat*p.csoil*crossect[j];
    }
    return fin;
  }

  pure double fluxin_sum(in ref Parameters p){
    double[] fin = fluxin(p);
    return sum(fin);
  }

  pure double fluxside_sum(){
    return sum(s.fluxside);
  }

  pure double fluxtop_a_sum(){
    return sum(s.fluxtop_a);
  }

  pure double fluxtop_d_sum(){
    return sum(s.fluxtop_d);
  }

  /**
   *  fills the 4 value array used for checking equilibrium
   */
  void fill_summary(ref double[4] s){
    s[0] = storage_sum();
    s[1] = fluxside_sum();
    s[2] = fluxtop_d_sum();
    s[3] = fluxtop_a_sum();
  }

  Matrix copy_matrix(string x){
    Matrix copy;
    copy.initialise(ny,nr);
    for(auto i = 0; i < ny; i++){
      for(auto j = 0; j < nr; j++){
	if(x == "cair"){
	  copy[i,j] = s.cair[i,j];
	}else if(x == "cwat"){
	  copy[i,j] = s.cwat[i,j];
	}else if(x == "nair"){
	  copy[i,j] = s.nair[i,j];
	}else if(x == "nwat"){
	  copy[i,j] = s.nwat[i,j];
	}else if(x == "towater"){
	  copy[i,j] = towater[i,j];
	}
      }
    }
    return copy;
  }
  
};

unittest
{
  Parameters p;
  p.set_defaults();
  Tree t = new Tree(p);
  assert(t.s.time < 0.1);
}
