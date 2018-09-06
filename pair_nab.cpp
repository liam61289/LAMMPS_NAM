/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_nab.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairNAB::PairNAB(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  nmax = 0;
  rho = NULL;
  fp = NULL;

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
  
  // -------- Modification CPR ------------
  if (force->newton_pair) {
    error->all(FLERR,"Set newton off to use pair_nab");
  }
  // -------- Modification CPR ------------
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairNAB::~PairNAB()
{
  memory->destroy(rho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete [] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
    }
    memory->sfree(funcfl);
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }

  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairNAB::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  // -------- Modification CPR ------------
  double H; // Hopping integral
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz; // Relative velocity
  double vdotr,invr2,invrhoi2,invrhoj2,invrhoi,invrhoj;
  double rhop,rhop2;
  double fx,fy,fz;
  double dostemp,hotemp;
  double gamma;
  // -------- Modification CPR ------------

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **f = atom->f;
  // -------- Modification CPR ------------
  double **v = atom->v;
  double **ffric = atom->ffric;
  double *ldamp = atom->ldamp;
  // -------- Modification CPR ------------
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    //printf ("atom: %i has %i neighbours (nlocal %i):", i, jnum, nlocal);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      //if (j < nlocal) {printf (" %i ," , j);}

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      
      //if (j == 904){printf ("%i %i: %f %f %f: %f %f %f\n",j,i,x[j][0],x[j][1],x[j][2],x[i][0],x[i][1],x[i][2]);}

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        // -------- Modification CPR ------------
        // Calculate sum sq hopping integrals for all atoms and
        // load into rho[]
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        //rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        H = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        
        rho[i] += H*H;
        /* if (i == 42) {
            printf ("H %f: %f\n", sqrt(rsq),H);
        } */
        //rho[i] += 1.0;
        //printf ("i %i %e:\n", i, rho[i]);
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          //rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          H = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          rho[j] += H*H;

          //rho[j] += 1.0;
          //printf ("j %i %e:\n", i, rho[i]);
        }
        // -------- Modification CPR ------------
      }
    }
  }
  //for (ii = 0; ii < inum; ii++) {
  //  i = ilist[ii];
  //  printf ("%i %i %f:\n",ii, i, rho[i]);
  //}

  // communicate and sum densities

  //if (newton_pair) comm->reverse_comm_pair(this);
  
  //comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  /*for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);*/

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    // -------- Modification CPR ------------
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    // -------- Modification CPR ------------
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // -------- Modification CPR ------------
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;      
      delvx = v[j][0] - vxtmp;
      delvy = v[j][1] - vytmp;
      delvz = v[j][2] - vztmp;
      // -------- Modification CPR ------------

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        
        // -------- Modification CPR ------------
        
        // First get gradient of hopping integral at 
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        rhop = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        rhop2 = rhop*rhop;
        /* if (i == 42) {
            printf ("H' %f: %f\n", sqrt(rsq),rhop);
        } */
        vdotr = delvx*delx + delvy*dely + delvz*delz;
        
        invr2 = 1.0/rsq;
        invrhoi2=1.0/rho[i];
        invrhoj2=1.0/rho[j];
        invrhoi=sqrt(invrhoi2);
        invrhoj=sqrt(invrhoj2);
        //if (rho[j] == 0){printf ("%i %i: %f %f %f: %f %f %f\n",j,i,x[j][0],x[j][1],x[j][2],x[i][0],x[i][1],x[i][2]);}
        fx = invrhoi*invrhoj*vdotr*invr2*delx*rhop2;
        fy = invrhoi*invrhoj*vdotr*invr2*dely*rhop2;
        fz = invrhoi*invrhoj*vdotr*invr2*delz*rhop2;
        //f[i][0] += fx;
        //f[i][1] += fy;
        //f[i][2] += fz;
        ffric[i][0] += fx;
        ffric[i][1] += fy;
        ffric[i][2] += fz;
        ldamp[i] += invrhoi*invrhoj*rhop2;
        //ffric[i][0] = rho[i];
        //ffric[i][1] += invrhoj;
        //ffric[i][2] += vdotr;
        //printf ("%i %e:", i, ffric[i][3]);
        if (newton_pair || j < nlocal) {
          fx = invrhoj*invrhoi*vdotr*invr2*delx*rhop2;
          fy = invrhoj*invrhoi*vdotr*invr2*dely*rhop2;
          fz = invrhoj*invrhoi*vdotr*invr2*delz*rhop2;
          //f[j][0] -= fx;
          //f[j][1] -= fy;
          //f[j][2] -= fz;
          ffric[j][0] -= fx;
          ffric[j][1] -= fy;
          ffric[j][2] -= fz;
          ldamp[j] += invrhoj*invrhoi*rhop2;
          //ffric[j][0] += invrhoj2;
          //ffric[j][1] += invrhoj;
          //ffric[j][2] += rhop2;
        }
        // -------- Modification CPR ------------

        /*coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        fpair = -psip*recip;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
        // -------- Modification CPR ------------
        ffric[i][0] += delx*fpair;
        ffric[i][1] += dely*fpair;
        ffric[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          ffric[j][0] -= delx*fpair;
          ffric[j][1] -= dely*fpair;
          ffric[j][2] -= delz*fpair;
        }
        // -------- Modification CPR ------------
        */

        if (eflag) evdwl = phi;
        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNAB::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNAB::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairNAB::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2],funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *)
      memory->srealloc(funcfl,nfuncfl*sizeof(Funcfl),"pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file,arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i = j

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (i == j) {
        setflag[i][i] = 1;
        map[i] = ifuncfl;
        //atom->set_mass(i,funcfl[ifuncfl].mass);
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNAB::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNAB::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++){
      	cutmax = MAX(cutmax,funcfl[m].cut);
      	a = MAX(a,funcfl[m].a);
       strcpy(crys,funcfl[m].crys);
    }
  } else if (setfl){
       cutmax = setfl->cut;
       a = setfl->a;
	strcpy(crys,setfl->crys);
   }
    else if (fs){
       cutmax = fs->cut;
       a = fs->a;
	strcpy(crys,fs->crys);
    }
  
  cutforcesq = cutmax*cutmax;
  
  force->D=calculate_D();

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairNAB::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl-1];

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = fopen(filename,"r");
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  int tmp;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %lg %s",&tmp,&file->mass,&file->a,file->crys);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->crys,4,MPI_CHAR,0,world);
  MPI_Bcast(&file->mass,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->a,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  memory->create(file->frho,(file->nrho+1),"pair:frho");
  memory->create(file->rhor,(file->nr+1),"pair:rhor");
  memory->create(file->zr,(file->nr+1),"pair:zr");

  if (me == 0) grab(fptr,file->nrho,&file->frho[1]);
  MPI_Bcast(&file->frho[1],file->nrho,MPI_DOUBLE,0,world);

  if (me == 0) grab(fptr,file->nr,&file->zr[1]);
  MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) grab(fptr,file->nr,&file->rhor[1]);
  MPI_Bcast(&file->rhor[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairNAB::file2array()
{
  int i,j,k,m,n;
  int ntypes = atom->ntypes;
  double sixth = 1.0/6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1) * file->dr);
    rhomax = MAX(rhomax,(file->nrho-1) * file->drho);
    mass = MAX(mass,file->mass);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r,p,cof1,cof2,cof3,cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      frho[n][m] = cof1*file->frho[k-1] + cof2*file->frho[k] +
        cof3*file->frho[k+1] + cof4*file->frho[k+2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      rhor[n][m] = cof1*file->rhor[k-1] + cof2*file->rhor[k] +
        cof3*file->rhor[k+1] + cof4*file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
        r = (m-1)*dr;

        p = r/ifile->dr + 1.0;
        k = static_cast<int> (p);
        k = MIN(k,ifile->nr-2);
        k = MAX(k,2);
        p -= k;
        p = MIN(p,2.0);
        cof1 = -sixth*p*(p-1.0)*(p-2.0);
        cof2 = 0.5*(p*p-1.0)*(p-2.0);
        cof3 = -0.5*p*(p+1.0)*(p-2.0);
        cof4 = sixth*p*(p*p-1.0);
        z2r[n][m] = cof1*ifile->zr[k-1] + cof2*ifile->zr[k] +
          cof3*ifile->zr[k+1] + cof4*ifile->zr[k+2];
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairNAB::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairNAB::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairNAB::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while (ptr = strtok(NULL," \t\n\r\f")) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

double PairNAB::single(int i, int j, int itype, int jtype,
                       double rsq, double factor_coul, double factor_lj,
                       double &fforce)
{
  int m;
  double r,p,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  double *coeff;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);

  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp[i]*rhojp + fp[j]*rhoip + phip;
  fforce = -psip*recip;

  return phi;
}

/* ---------------------------------------------------------------------- */

int PairNAB::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairNAB::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) rho[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairNAB::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairNAB::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairNAB::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   swap fp array with one passed in by caller
------------------------------------------------------------------------- */

void PairNAB::swap_eam(double *fp_caller, double **fp_caller_hold)
{
  double *tmp = fp;
  fp = fp_caller;
  *fp_caller_hold = tmp;
}

double PairNAB::calculate_D()
{ //currently only works for materials of one atom type

  double D;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double H; // Hopping integral
  double p, rsq;
  int m, itype, jtype;
  int *type = atom->type;
  itype = type[0];
  jtype = type[0];

  	//calculate D, needed by fix_nab - modification Liam
	// Set up distance to each group of nearest neighbours.
	// Need this for each crystal structure
	int size =0;
	double *rd;
	double *rhod;
	double *num;
	
	if(strcmp(crys,"FCC")==0){
		size = 3;

		rd =  new double[size];
		rd[0]=(0.5*sqrt(2)*a); 
		rd[1]=(a);
		rd[2]=0.5*(sqrt(6)*a);

		rhod =  new double[size];		
		rhod[0]=0;
		rhod[1]=0;
		rhod[2]=0;
		
		num =  new double[size];
		num[0] = 12; 
		num[1] = 6;
		num[2] = 24;
	}
	else if(strcmp(crys,"BCC")==0){
		size = 3;

		rd =  new double[size];
		rd[0]=(0.5*sqrt(2)*a); 
		rd[1]=(a);
		rd[2]=0.5*(sqrt(6)*a);

		rhod =  new double[size];		
		rhod[0]=0;
		rhod[1]=0;
		rhod[2]=0;
		
		num =  new double[size];
		num[0] = 8; 
		num[1] = 6;
		num[2] = 24;

	}
	else if(strcmp(crys,"SC")==0){
		size = 3;

		rd =  new double[size];
		rd[0]=(a); 
		rd[1]=(sqrt(2)*a);
		rd[2]=(2*a);

		rhod =  new double[size];		
		rhod[0]=0;
		rhod[1]=0;
		rhod[2]=0;
		
		num =  new double[size];
		num[0] = 6; 
		num[1] = 12;
		num[2] = 6;
	}
	else if(strcmp(crys,"HCP")==0){
		error->all(FLERR,"HCP crystal structure not currently accepted");

	}
	else {
	    error->all(FLERR,"Unknown crystal structure");
	}

  	// loop to calculate rho for set up crystal structure
	// this loop is generic i.e. works for all specified crystals
   	for (int jj = 0; jj < size; jj++) {
     		rsq = rd[jj]*rd[jj];
             	if (rsq < cutforcesq) {
       		p = sqrt(rsq)*rdr + 1.0;
        		m = static_cast<int> (p);
        		m = MIN(m,nr-1);
        		p -= m;
        		p = MIN(p,1.0);
        		// -------- Modification CPR ------------
        		// Calculate sum sq hopping integrals for all atoms and
       		// load into rho[]
        		coeff = rhor_spline[type2rhor[jtype][itype]][m];
        		H = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        		rhod[jj] += H*H;
        		// -------- Modification CPR ------------
      		}
    	}
	// loop to add up rho for all nearest neigbours
	double Dtemp=0;
   	for (int jj = 0; jj < size; jj++) {
		Dtemp += num[jj]*rhod[jj];
	}
	D = 1/sqrt(Dtemp);

return D;
}