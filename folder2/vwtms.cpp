#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <fstream>

#include "../space/Space.h"
#include "../space/globle.h"
#include "../space/grid_space.h"
#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
//void to take a van der Waals epsmap and expand it into a
//molecular volume map..
//procedure is to first generate a list of accessible points
//whose collective volume may contain any midpoint in the box.
//next s}// these points to be mapped into indexing arrays...
//{ take the vanderwaals epsmap and grow it out to the molecular
//epsmap, checking against the set of accessible points..
//
//can't do the vw surface by projecting midpoints to accessible surface
//and testing if that point is in or outside the accessible volume
//(20 Oct 92)
//(but can do the contact part of MS that way, Nov 93)
//
//program to expand the vw surface into Richards' molecular surface
//(S. Sridharan May 1994)
//
//2011-05-17 Parameters transferred via modules qlog and pointers
//
//int cbn1(1),cbn2(1),cbal(1) !2011-05-16 Arrays declared
//int iab1(1),iab2(1),icume(1) !in pointers module
//pointer (i_egrid,egrid) 2011-05-16 transfered to pointers.f95
//#####################################################################

void vwtms;

bool outcb[-2:2,-2:2,-2:2];
character(80) :: line;
int nbra[1000];
//2011-05-17 bndeps is declared in pointers module
int eps[6],nt;
int itmp[6],dim1,cont;
int iaprec,objecttype,imap[4,6],dim,sign,epsdim,isign;
int kind,eps2[6],j3;
bool remov;
character(96) :: strtmp;
SGrid <float> itemp, rtemp,xq,dxyz,dr123,dx123,u123;
SGrid <float> goff[6],xg,x1xyz,s123,xxyyzz;
SGrid <int> ixyz,it123,ixyz2,jxyz;
type(int_coord), allocatable :: ibgrd_temp(:);
bool out,nbe[0:6],intb,exists,flag;
//2011-05-27 Declarations added due to IMPLICIT NONE
float offf,cbln,cba,del,dis,dmn,dist,ds2,off,dsr,r0a;
float s1,s2,s3,x1,radp,prbrd12,prbrd22,rsm,rsm1,prbrd2;
int ia,i,iac,ibgp,ii,iext,iarv,iacv,iiii,imedia;
int iiord,ios,it1,it2,it3,iv,ix2,iy2,iz2,ix,iy,iz;
int iii,iord,j,jj,jjj,jx,jy,jz,limu,liml,kk,m,k,n;
int mr,mpr,n2,nacc,natm,nmmt,ndv,ncms,ncav,nnn,nnbr;
int nn,ntt,nmt,n1,Nar,iqqq;

//2011-05-17 Arrays allocated by ordinary F95 allocate statement
allocate(ibnd(ibmx),r0(iNatom),r02(iNatom),rs2(iNatom));
allocate(ast(iNatom));
allocate(bndeps[iGrid][iGrid][iGrid][2]);

offf=(iGrid+1.)/2.;
epsdim=iNatom+iNObject+2;

//imap maps from midpoint position to iepsmap entry positions
//2011-05-17 Changed to array operations
imap=0;

imap[1][4]=-1;
imap[2][5]=-1;
imap[3][6]=-1;
imap[4][1]=1;
imap[4][2]=2;
imap[4][3]=3;
imap[4][4]=1;
imap[4][5]=2;
imap[4][6]=3;

//2011-05-17 Changed to array operations
outcb=true;
outcb[-1:1][-1:1][-1:1]=false;
nbe=false;
nbe[1:5]=true;

goff=coord(0.,0.,0.);
off=0.5/fScale;

//2011-05-17 Changed to coord derived variable type
goff[1].nX=off;
goff[2].nY=off;
goff[3].nZ=off;
goff[4].nX=-off;
goff[5].nY=-off;
goff[6].nZ=-off;

radpmax=max(fRadPrb(1),fRadPrb(2)) //transferred via qlog module;

//convertion from grid to float coordinates(can also use routine
//gtoc)
//2011-05-17 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
x1=1.0/fScale;
x1xyz=cOldMid-(0.5*x1*(iGrid+1));

//find extrema

//find global extrema
cMin=coord(6000,6000,6000);
cMax=coord(-6000,-6000,-6000);
for(ii=1;ii<=iNObject;ii++){
cMin=min(cMin][sLimObject[ii]%min);
cMax=max(cMax][sLimObject[ii]%max);
}// do

//find vanderwaals boundary
n=0;
nn=0;
nmt=0;
nmmt=0;

//NB change limits to those of the molecule.
//set for iEpsMap NOT equal to unity
for(k=LimEps%min.nZ+1;k<=LimEps%max.nZ-1;k++){
for(j=LimEps%min.nY+1;j<=LimEps%max.nY-1;j++){
for(i=LimEps%min.nX+1;i<=LimEps%max.nX-1;i++){
//one distinguishes between internal,external,
//internal bgp and external bgp
iext=0;
ibgp=0;

//2011-05-17 Changed to iEpsMap to int_coord derived type
itmp[1]=iabs(iEpsMap[i][j][k].nX)/epsdim;
itmp[2]=iabs(iEpsMap[i][j][k].nY)/epsdim;
itmp[3]=iabs(iEpsMap[i][j][k].nZ)/epsdim;
itmp[4]=iabs(iEpsMap[i-1][j][k].nX)/epsdim;
itmp[5]=iabs(iEpsMap[i][j-1][k].nY)/epsdim;
itmp[6]=iabs(iEpsMap[i][j][k-1].nZ)/epsdim;

if(itmp[1]==0) iext=1;
if(itmp[1]!=itmp[6]) ibgp=1;

for(cont=2;cont<=6;cont++){
if(itmp[cont]==0) iext=1;
if(itmp[cont]!=itmp[cont-1]) ibgp=1;
}// do

//assignement of right values to bndeps according to the
//point nature
//from now iBoundNum is the total number of internal and
//external boundary grid points
if (ibgp>0) {
n=n+1;
bndeps[i][j][k][1]=n;
bndeps[i][j][k][2]=iext;
if (iext>0) nn=nn+1;
ibnd(n)=int_coord(i,j,k);
}
else{
bndeps[i][j][k][1]=0;
bndeps[i][j][k][2]=0;
}// if

if (bDebug) {
nt=0;
//passed from mod to dim, I need the medium
if((iEpsMap[i][j][k].nX/epsdim)>0)nt=nt+1;
if((iEpsMap[i][j][k].nY/epsdim)>0)nt=nt+1;
if((iEpsMap[i][j][k].nZ/epsdim)>0)nt=nt+1;
if((iEpsMap[i-1][j][k].nX/epsdim)>0)nt=nt+1;
if((iEpsMap[i][j-1][k].nY/epsdim)>0)nt=nt+1;
if((iEpsMap[i][j][k-1].nZ/epsdim)>0)nt=nt+1;

if (nbe[nt].neqv.((iext==1)&&(ibgp==1))){
cout <<"PROBLEMS1 " << i << j << k << endl;

//2011-05-17 converts grid to float coordinates
//using operations on coord and int_coord type
//variables defined in module
//operators_on_coordinates
itemp=Int2Float(int_coord(i,j,k));
rtemp=((itemp-offf)/fScale)+cOldMid;

cout <<rtemp << endl;
cout <<iEpsMap[i << j << k] << iEpsMap[i-1 << j << k].nX << iEpsMap[i << j-1 << k].nY << iEpsMap[i << j << k-1].nZ << endl;
}// if
}// if
}// do
}// do
}// do

iBoundNum=n;
iBoundNumsurf=nn;
nn=0;

if (bVerbose) {
cout <<"boundary points facing continuum solvent= " << iBoundNumsurf << endl;
cout <<"total number of boundary points before elab.= " << iBoundNum << endl;
}// if

if (iBoundNum>ibmx) {
cout <<"iBoundNum= " << iBoundNum << " is greater than ibmx = " << ibmx << endl;
cout <<"increase ibmx in vwtms.f" << endl;
stop;
}// if

if (bDebug) {
open(52,file='bgpPrima.log',form='formatted');

//grid coordinates
write(52,*)iGrid,iBoundNum;
for(iiii=1;iiii<=iBoundNum;iiii++){
write(52,*)iiii,ibnd(iiii);
}// do
close (52);
}// if //end bDebug;

//2011-05-17 Arrays allocated by ordinary F95 allocate statement
if (radpmax<1.e-6) {
allocate(ibgrd(iBoundNum));

//2011-05-17 Changed to array operations, but keeping some
//assignment in a cycle due to array size mismatch
ast=0;
for(i=1;i<=iBoundNum;i++){
ibgrd(i)=ibnd(i);
}// do
}
else{
if (!bOnlyMol&&fRadPrb(1)!=fRadPrb(2)) {
for(i=1;i<=iNatom;i++){
//2011-05-17 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
xq=xn1(i);
r0a=sDelPhiPDB[i].radius+fRadPrb(1);

for(ii=1;ii<=iNObject;ii++){
strtmp=dataobject[ii][1];
read(strtmp(16:18),*)kind;
if (strtmp(1:4)!='is a'&&kind!=2) {
if (((xq-sDelPhiPDB[i].radius).vandlt.&sLimObject[ii]%max)&&((xq+&sDelPhiPDB[i].radius).vandgt.sLimObject[ii]%min)) {

call distobj(xq,dist,dxyz,ii,0.,true);

//only full buried atoms have the proberadius
//changed to fRadPrb(2)
if (dist<-sDelPhiPDB[i].radius) r0a=sDelPhiPDB[i].radius+fRadPrb(2);
}// if
}// if
}// do

r0(i)=r0a;
r02(i)=r0a*r0a;
rs2(i)=0.99999*r02(i);
}// do
}
else{
for(i=1;i<=iNatom;i++){
r0a=sDelPhiPDB[i].radius+fRadPrb(1);
r0(i)=r0a;
r02(i)=r0a*r0a;
rs2(i)=0.99999*r02(i);
}// do
}// if

if(iacs) write(6,'(a21)') " opening surface file";
if(iacs) open(40,file='hsurf2.dat');

//make a list of accessible points..,expos. all scaling of grid
//points will be done to thses points..
prbrd12=fRadPrb(1)*fRadPrb(1);
prbrd22=fRadPrb(2)*fRadPrb(2);

//calculate an identity for this conformation
rsm=0.0;
if (bOnlyMol) {
for(i=1;i<=iNatom;i++){
rsm=rsm+sDelPhiPDB[i].radius*sum(abs(xn1(i)));
}// do
}// if

inquire(file='ARCDAT',exist=exists);
flag=true;
if (exists&&bOnlyMol) {
open(1,file='ARCDAT',form='unformatted',status='old',iostat=ios);
read(1)natm,radp,nacc,rsm1;

//maybe it could be improved taking into account objects
//(Walter 1999)
if (natm==iNatom&&radp==fRadPrb(1)&& abs(rsm1-rsm)<1.0e-8) {
if(bVerbose)cout <<"reading accessible surface arcs data from file ARCDAT" << endl;
extot=nacc;
allocate(expos(extot));

//2011-05-19 Subroutine arcio is too short and called
//only from this void, thus no meaning to be
//separate piece of code
read(1)ast;
read(1) expos;
if(bVerbose) cout <<"no. of arc points read = " << nacc << endl;
close (1);
flag=false;
}
else{
close (1);
}// if
}// if

if (flag) {
call sas(xn1,extot);

if (extot>0&&bOnlyMol) {
//maybe can be improved by taking into account objects
cout <<"writing accessible surface arcs data to ARCDAT" << endl;
open(1,file='ARCDAT',form='unformatted',iostat=ios);

if (ios!=0) {
cout << "error opening ARCDAT file" << endl;
stop;
}// if

write(1)iNatom,fRadPrb(1),extot,rsm;

if (iNatom==0) {
write(1)0;
write(1) expos;
}
else{
write(1)ast;
write(1)expos;
}// if

close(1);
}// if

if (bDebug) {
open(52,file='Vertices.txt',form='formatted');
for(iiii=1;iiii<=extot;iiii++){
write (52,*) expos(iiii);
}// do
close (52);
}// if //end bDebug;

}// if

del=1./fScale;
del=max(del,radpmax);
cbln=fRMax+del;

call cubedata(2.0,cbln);

dim=(lcb+1)*(mcb+1)*(ncb+1);

//2011-05-17 Array allocation is done by ordinary F95 ALLOCATE
//statement Array allocated as 3D as in the cube void
allocate(cbn1(dim),cbn2(dim));

dim1=27;
if ((iNObject-numbmol)>0) dim1=max(dim,27);
allocate(cbal(dim1*(iNatom+iNObject-numbmol)));

call cube(xn1,fRadPrb(1),cbn1,cbn2);

//link the accessible points into iab1 and iab2
call indverdata(radpmax,fScale);

cba=1./grdi;
NNN=(lcb1+1)*(mcb1+1)*(ncb1+1);
allocate(iab1(NNN),iab2(NNN));
allocate(icume(extot));

if(bVerbose)write(6,"(a,f8.3)")' grid for indexing accessible points = ',cba;

call indver(extot,iab1,iab2);

//write out surface data
if (iacs) {
line= ' ';
line(1:6)='ATOM ';
line(14:14)='O';
line(18:20)='SP ';

for(i=1;i<=extot;i++){
xg=expos(i);
iv=1;

//2011-05-24 Subroutine watput is small and called only
//once, thus the code transfered into calling void
write(line(7:11),'(i5)')i;
write(line(24:26),'(i3)')iv;

//2011-05-24 In original coordinates from array xo were
//written into th line.
//Logic of this piece, however, requires coordinates
//from xg to be written
write(line(31:54),'(3(f8.3))')xg;
write(40,'(a80)') line;
}// do
close (40);
}// if

//now start the expansion
//m1= the number of boundary points removed from list
ncav=0;
n1=1;
n2=iBoundNum;

//m= number of new boundary elements..
mpr=100000;
ndv=0;

D100: do
m=0;
mr=0;

for(i=n1;i<=n2;i++){
ixyz=ibnd(i);
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;

//considering both internal and external b.g.p.
if (bndeps[ix][iy][iz][1]!=0) {

//still has to be considered what is external and
//what internal!!!!!WWW
//remov is true if it is an internal midpoint close
//to an interface where a molecule is present
//(expansion has to take place also in objects)
remov=false;

//tengo il mod perche' deve pr}//ere solo punti in
//atomi
eps[1]=mod(iEpsMap[ix][iy][iz].nX][epsdim);
eps[2]=mod(iEpsMap[ix][iy][iz].nY][epsdim);
eps[3]=mod(iEpsMap[ix][iy][iz].nZ][epsdim);
eps[4]=mod(iEpsMap[ix-1][iy][iz].nX][epsdim);
eps[5]=mod(iEpsMap[ix][iy-1][iz].nY][epsdim);
eps[6]=mod(iEpsMap[ix][iy][iz-1].nZ][epsdim);

remov=((eps[1]>1&&eps[1]<=iNatom+1)||(eps[2]>1&&eps[2]<=iNatom+1));
remov=((eps[3]>1&&eps[3]<=iNatom+1)||(eps[4]>1&&eps[4]<=iNatom+1))||remov;
remov=((eps[5]>1&&eps[5]<=iNatom+1)||(eps[6]>1&&eps[6]<=iNatom+1))||remov;

//da farsi solo se pores eps2 contiene il mezzo
eps2[1]=(iEpsMap[ix][iy][iz].nX/epsdim);
eps2[2]=(iEpsMap[ix][iy][iz].nY/epsdim);
eps2[3]=(iEpsMap[ix][iy][iz].nZ/epsdim);
eps2[4]=(iEpsMap[ix-1][iy][iz].nX/epsdim);
eps2[5]=(iEpsMap[ix][iy-1][iz].nY/epsdim);
eps2[6]=(iEpsMap[ix][iy][iz-1].nZ/epsdim);

//cWWW there is still an issue in case there are both
//molecules and objects: since parent object of
//reentrant points is only known in sclbp, filling
//reentrant regions due to molecules in objects might
//fail
remov=remov&&(numbmol>0);

//2011-05-17 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
xg=(x1*Int2Float(ixyz))+x1xyz;

D200: do j=1,6;
//essere in poro ==> eps2=0 and eps >0
if (eps[j]==0||(remov&&eps[j]>iNatom+1)||(eps2[j]==0&&eps[j]>0)) {
prbrd2=prbrd22;
if (eps[j]==0||eps2[j]==0) prbrd2=prbrd12;

//add midpoint offset to grid point..
s123=xg+goff[j];

//determine if this virgin midpoint is in or out
//2011-05-18 mn(x,y,z) and grdi were
//assigned values in INDVER void
//now coord type variable mnxyz is declared
//in pointers module and float grdi declared
//and thus accessible in qlog module
xxyyzz=(s123-mnxyz)*grdi;
jxyz=Float2Int(xxyyzz);
jx=jxyz.nX;
jy=jxyz.nY;
jz=jxyz.nZ;

//2011-05-18 Indexes lcb1, mcb1, ncb1 are
//transfred via qlog module and are set in
//INDVER void
if ((jxyz.vorle.0)||(jxyz.vorge.lmncb1)) {
cout <<"midpoint out of cube" << endl;
write(6,'(2i5,3f8.3,3i6,3i8)')i,j,xxyyzz,jxyz,lmncb1;
cout <<iEpsMap[ix << iy << iz].nX << endl;
cout <<iEpsMap[ix << iy << iz].nY << endl;
cout <<iEpsMap[ix << iy << iz].nZ << endl;
cout <<iEpsMap[ix-1 << iy << iz].nX << endl;
cout <<iEpsMap[ix << iy-1 << iz].nY << endl;
cout <<iEpsMap[ix << iy << iz-1].nZ << endl;
}// if

dmn=1000.;
iacv=0;

//2011-05-18 Repeating piece of the code now
//is in the separate file
include 'inc_epsmakmod_vwtms.inc';

//-1,0,0
jx=jx-1;

include 'inc_epsmakmod_vwtms.inc';

//1,0,0
jx=jx+2;

include 'inc_epsmakmod_vwtms.inc';

//0,-1,0
jx=jx-1;
jy=jy-1;

include 'inc_epsmakmod_vwtms.inc';

//0,1,0
jy=jy+2;

include 'inc_epsmakmod_vwtms.inc';

//0,0,-1
jy=jy-1;
jz=jz-1;

include 'inc_epsmakmod_vwtms.inc';

//0,0,1
jz=jz+2;

include 'inc_epsmakmod_vwtms.inc';

//nn=2
//1,0,1
jx=jx+1;

include 'inc_epsmakmod_vwtms.inc';

//-1,0,1
jx=jx-2;

include 'inc_epsmakmod_vwtms.inc';

//0,1,1
jx=jx+1;
jy=jy+1;

include 'inc_epsmakmod_vwtms.inc';

//0,-1,1
jy=jy-2;

include 'inc_epsmakmod_vwtms.inc';

//-1,-1,0
jz=jz-1;
jx=jx-1;

include 'inc_epsmakmod_vwtms.inc';

//1,-1,0
jx=jx+2;

include 'inc_epsmakmod_vwtms.inc';

//1,1,0
jy=jy+2;

include 'inc_epsmakmod_vwtms.inc';

//-1,1,0
jx=jx-2;

include 'inc_epsmakmod_vwtms.inc';

//-1,0,-1
jz=jz-1;
jy=jy-1;

include 'inc_epsmakmod_vwtms.inc';

//1,0,-1
jx=jx+2;

include 'inc_epsmakmod_vwtms.inc';

//0,1,-1
jx=jx-1;
jy=jy+1;

include 'inc_epsmakmod_vwtms.inc';

//0,-1,-1
jy=jy-2;

include 'inc_epsmakmod_vwtms.inc';

//nn=3
//-1,-1,-1
jx=jx-1;

include 'inc_epsmakmod_vwtms.inc';

//1,-1,-1
jx=jx+2;

include 'inc_epsmakmod_vwtms.inc';

//1,1,-1
jy=jy+2;

include 'inc_epsmakmod_vwtms.inc';

//-1,1,-1
jx=jx-2;

include 'inc_epsmakmod_vwtms.inc';

//-1,1,1
jz=jz+2;

include 'inc_epsmakmod_vwtms.inc';

//1,1,1
jx=jx+2;

include 'inc_epsmakmod_vwtms.inc';

//1,-1,1
jy=jy-2;

include 'inc_epsmakmod_vwtms.inc';

//-1,-1,1
jx=jx-2;

include 'inc_epsmakmod_vwtms.inc';

//it might be in the contact region
find;
//the closest atom surface
//2011-05-18 Using operations on coord and
//int_coord type variables defined
//in module operators_on_coordinates
it123=Float2Int((s123-xyzo)*cbai);

dmn=100.;
iac=0;
nnbr=0;
lmncb=int_coord(lcb,mcb,ncb);

if((it123.vorlt.0)||(it123.vorgt.lmncb)){
//if the bgp is outside the cube,
//probably it is due to some object
for(ii=iNObject;ii<=1,-1;ii++){
strtmp=dataobject[ii][1];
read(strtmp(16:18),*)kind;
if(strtmp(1:4)!='is a'&&kind!=2) {
if ((s123.vandle.(sLimObject[ii]%max+&x1))&&(s123.vandgt.(sLimObject[ii]%min-x1))) {
nnbr=nnbr+1;
nbra[nnbr]=ii+iNatom;
liml=0;
limu=0;
}// if
}// if
}// do

if(liml!=0||limu!=0) cout <<"a bgp close to nothing" << endl;
}
else{
//2011-05-19 Changed 1d array to 3d array as
//in cube void
liml=cbn1(it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ);
limu=cbn2(it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ);
}// if

iaprec=0;
DOKK: do kk=liml,limu;
ia=cbal(kk);
if (ia==0) cout <<"problems with cube" << endl;

if (ia<=iNatom&&ia>0) {
if (ast(ia)==0) {
nnbr=nnbr+1;
nbra[nnbr]=ia;
}// if
}
else{
if (ia!=iaprec&&eps[j]==0) {
//different from sclbp, I have to
//consider the object only if the
//midpoint is external O SE C'E' UN
//PORO E VEDERE IN SEGUITO
iaprec=ia;

//assuming any object is not buried
nnbr=nnbr+1;
nbra[nnbr]=ia;
}// if
}// if
}// do DOKK;

DOII: do ii=1,nnbr;
ia=nbra[ii];

if (ia>iNatom) {
iii=ia-iNatom;
//2011-05-18 Using operations on coord
//and int_coord type variables defined
//in module operators_on_coordinates
xq=s123;

//check if bgp is inside VdW of object
call distobj(xq,dist,dxyz,iii,0.0,false);

//an object can compete with an atom for
//a midpoint only if this midpoint
//is out of the object itself
if (dist>=0.&&dist<dmn) {
dmn=dist;
iac=ia;
dr123=(-dxyz)*(fRadPrb(1)-dist);
}// if
}
else{
dx123=s123-xn1(ia);
ds2=dx123.dot.dx123;
dis=sqrt(ds2)-sDelPhiPDB[ia].radius;

//dis= distance to atom surface
if (dis<dmn){
dmn=dis;
iac=ia;
}// if
}// if
}// do DOII;

if (iac==0) {
if (bDebug) {
cout <<"bgp:" << i << " might be a cavity point" << ix << iy << iz << endl;
cout <<"midpoint" << j << " in position [A]" << s1 << s2 << s3 << endl;
cout <<"it1:" << it1 << " it2:" << it2 << " it3:" << it3 << endl;
}// if

ncav=ncav+1;

//possibly a cavity point
}
else{
//check to see if it is in the contact
//region of that atom or object by
//projecting it to the atom's acc surface
//and checking against the acc volumes of
//nearby atoms
if (iac<=iNatom) {
dr123=s123-xn1(iac);
dsr=sqrt(dr123.dot.dr123);
u123=xn1(iac)+((r0(iac)*dr123)/dsr);
}
else{
u123=s123-dr123;
}// if

it123=Float2Int((u123-xyzo)*cbai);

liml=cbn1(it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ);
limu=cbn2(it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ);

DLIM: do kk=liml,limu;
ia=cbal(kk);
if (ia<=iNatom) {
dx123=u123-xn1(ia);
ds2=dx123.dot.dx123;
flag=true;

if (ds2<rs2(ia)) {
flag=false;
exit DLIM;
}// if
}
else{
if (ia!=iac&&eps[j]==0) {
xq=u123;

call distobj(xq,dist,dxyz,ia-iNatom,fRadPrb(1),true);

flag=true;
if (dist<-1.e-6) {
flag=false;
exit DLIM;
}// if

//oriented distance from ext}//ed
//object surface if negative =>
//reentrant region
}// if
}// if
}// do DLIM;

//it is in the contact region. flag the
//midpoint so it is not checked again
//iac is atom number...NOT increased by 1

//2011-05-18 To get rid of goto 201
//statements above
if (flag) {
eps[j]=-iac;
eps2[j]=-iac;
cycle D200;
}// if
}// if

eps[j]=1 //eps = 1 means cavity or reentrant;

//remap iEpsMap
if (iac==0) {
//this is an assumption, still to deeply
//understand meaning of cavity here and to
//improve this choice!WWW
if (ia>0) {
imedia=iAtomMed[ia];
}
else{
cout <<"assigning arbitrary epsilon in cavity" << endl;
imedia=iAtomMed[1];
}// if
}
else{
imedia=iAtomMed[iac];
}// if

switch (imap[4][j]){
case 1;
iEpsMap[ix+imap[1][j]][iy+imap[2][j]][iz+imap[3][j]].nX=eps[j]+imedia*epsdim;
case 2;
iEpsMap[ix+imap[1][j]][iy+imap[2][j]][iz+imap[3][j]].nY=eps[j]+imedia*epsdim;
case 3;
iEpsMap[ix+imap[1][j]][iy+imap[2][j]][iz+imap[3][j]].nZ=eps[j]+imedia*epsdim;
case default;
cout <<"?????" << endl;
}// select;

eps2[j]=imedia;

//not assigning the owner but only the medium,
//the former will be assigned in the fScale
//routine

//check to see if the nearest neighbour status
//has been changed..
ix2=ix;
iy2=iy;
iz2=iz;

//if the nearest neighbour is a box boundary
//point { skip this since box boundary
//points can not also be dielctric boundary
//points
//2011-05-18 Multiple IFs replaced by SELECT
//CASE
switch (j){
case 1;
ix2=ix+1;
if(ix2==iGrid) cycle D200;
case 2;
iy2=iy+1;
if(iy2==iGrid) cycle D200;
case 3;
iz2=iz+1;
if(iz2==iGrid) cycle D200;
case 4;
ix2=ix-1;
if(ix2==1) cycle D200;
case 5;
iy2=iy-1;
if(iy2==1) cycle D200;
case 6;
iz2=iz-1;
if(iz2==1) cycle D200;
}// select;

//once again one distinguishes between
//internal,external,internal bgp and external
//bgp
iext=0;
ibgp=0;

//2011-05-18 Changed to i nt_coord derived type
itmp[1]=iabs(iEpsMap[ix2][iy2][iz2].nX)/epsdim;
itmp[2]=iabs(iEpsMap[ix2][iy2][iz2].nY)/epsdim;
itmp[3]=iabs(iEpsMap[ix2][iy2][iz2].nZ)/epsdim;
itmp[4]=iabs(iEpsMap[ix2-1][iy2][iz2].nX)/epsdim;
itmp[5]=iabs(iEpsMap[ix2][iy2-1][iz2].nY)/epsdim;
itmp[6]=iabs(iEpsMap[ix2][iy2][iz2-1].nZ)/epsdim;

if(itmp[1]==0) iext=1;
if(itmp[1]!=itmp[6]) ibgp=1;

for(cont=2;cont<=6;cont++){
if(itmp[cont]==0) iext=1;
if(itmp[cont]!=itmp[cont-1]) ibgp=1;
}// do

if (bDebug) {
nt=0;
if((iEpsMap[ix2][iy2][iz2].nX/epsdim)>0) nt=nt+1;
if((iEpsMap[ix2][iy2][iz2].nY/epsdim)>0) nt=nt+1;
if((iEpsMap[ix2][iy2][iz2].nZ/epsdim)>0) nt=nt+1;
if((iEpsMap[ix2-1][iy2][iz2].nX/epsdim)>0) nt=nt+1;
if((iEpsMap[ix2][iy2-1][iz2].nY/epsdim)>0) nt=nt+1;
if((iEpsMap[ix2][iy2][iz2-1].nZ/epsdim)>0) nt=nt+1;
if(nbe[nt].neqv.(ibgp==1&&iext==1)) {
cout <<"PROBLEMS3" << ix2 << iy2 << iz2 << endl;
}// if
}// if //end+bDebugging;

if ((ibgp==0)&&(bndeps[ix2][iy2][iz2][1]!=0)) {
//reset bndeps for that point (i.e. remove
//bgp flag).
//a bgp become internal
iBoundNumsurf=iBoundNumsurf-bndeps[ix2][iy2][iz2][2];
bndeps[ix2][iy2][iz2][1]=0;
bndeps[ix2][iy2][iz2][2]=0;
mr=mr+1;
}
else{
if (ibgp==1&&iext==0&&bndeps[ix2][iy2][iz2][2]==1){
//an ext bgp is turned into an internal
//bgp
iBoundNumsurf=iBoundNumsurf-1;
bndeps[ix2][iy2][iz2][2]=0;
}// if
}// if

if (ibgp==1&&bndeps[ix2][iy2][iz2][1]==0) {
//create a new boundary point..
m=m+1;
bndeps[ix2][iy2][iz2][1]=n2+m;
ibnd(n2+m)=int_coord(ix2,iy2,iz2);
bndeps[ix2][iy2][iz2][2]=iext;
iBoundNumsurf=iBoundNumsurf+bndeps[ix2][iy2][iz2][2];
}// if
}// if

//now jump to the next midpoint of the same grid
//point
}// do D200;

//remap iEpsMap in case there have been changes..
//(that is some 0's became -1's)
//in other words: midpoint must remain external to
//objects
for(jj=1;jj<=6;jj++){
//in this way I can deal with eps[jj]<0
isign=1;

//iord=owner of the midpoint jj before change or
//after eps=1
iiord=comp(iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]]][imap[4][jj]);
iord=mod(iiord,epsdim);

//the last changed for sure has not iord<0
//there can be iord<0 due to nearest neighbors
//already changed
if (iord<0) cycle;

//if it has changed at previous step, dont change
//anymore
if (eps[jj]<0) {
isign=-1;
if (iord==0) iord=1;
}// if

jjj=iabs(iiord)/epsdim;

switch (imap[4][jj]){
case 1;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nX=isign*(iord+jjj*epsdim);
case 2;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nY=isign*(iord+jjj*epsdim);
case 3;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nZ=isign*(iord+jjj*epsdim);
case default;
cout <<"?????" << endl;
}// select;

//left iord definition with mod since if it is <>
//0, it keeps its identity
}// do

//at this point one still can trace what changes have
//been made check to see if this is still a boundary
//point
//once again one distinguishes between
//internal,external,internal bgp and external bgp
iext=0;
ibgp=0;

itmp[1]=iabs(iEpsMap[ix][iy][iz].nX)/epsdim;
itmp[2]=iabs(iEpsMap[ix][iy][iz].nY)/epsdim;
itmp[3]=iabs(iEpsMap[ix][iy][iz].nZ)/epsdim;
itmp[4]=iabs(iEpsMap[ix-1][iy][iz].nX)/epsdim;
itmp[5]=iabs(iEpsMap[ix][iy-1][iz].nY)/epsdim;
itmp[6]=iabs(iEpsMap[ix][iy][iz-1].nZ)/epsdim;

if(itmp[1]==0) iext=1;
if(itmp[1]!=itmp[6]) ibgp=1;

for(cont=2;cont<=6;cont++){
if(itmp[cont]==0) iext=1;
if(itmp[cont]!=itmp[cont-1]) ibgp=1;
}// do

if (bDebug) {
nt=0;
if((iabs(iEpsMap[ix][iy][iz].nX)/epsdim)>0)nt=nt+1;
if((iabs(iEpsMap[ix][iy][iz].nY)/epsdim)>0)nt=nt+1;
if((iabs(iEpsMap[ix][iy][iz].nZ)/epsdim)>0)nt=nt+1;
if((iabs(iEpsMap[ix-1][iy][iz].nX)/epsdim)>0) nt=nt+1;
if((iabs(iEpsMap[ix][iy-1][iz].nY)/epsdim)>0) nt=nt+1;
if((iabs(iEpsMap[ix][iy][iz-1].nZ)/epsdim)>0) nt=nt+1;

if (nbe[nt].neqv.(ibgp==1&&iext==1)) {
cout <<"PROBLEMS4" << ix << iy << iz << endl;
cout <<"epsdim=" << epsdim << "ibgp=" << ibgp << "iext=" << iext << endl;
cout <<"itmp" << itmp << endl;
cout <<iEpsMap[ix << iy << iz].nX << endl;
cout <<iEpsMap[ix << iy << iz].nY << endl;
cout <<iEpsMap[ix << iy << iz].nZ << endl;
cout <<iEpsMap[ix-1 << iy << iz].nX << endl;
cout <<iEpsMap[ix << iy-1 << iz].nY << endl;
cout <<iEpsMap[ix << iy << iz-1].nZ << endl;
}// if
}// if

//if not now a boundary element change bndeps
if ((iext==0)||(ibgp==0)) {
iBoundNumsurf=iBoundNumsurf-bndeps[ix][iy][iz][2];
if(ibgp==1) bndeps[ix][iy][iz][2]=iext;

if (ibgp==0) {
bndeps[ix][iy][iz][1]=0;
bndeps[ix][iy][iz][2]=0;
mr=mr+1;
if(iext==1)cout <<"//!!born a new external point!!!" << endl;
}// if
}// if

}// if //if end for whether bndeps is nonzero;
}// do //next boundary point FINISH;

n1=n2+1;
n2=n2+m;
if(bVerbose) cout <<"bgp added m=" << m << " bgp removed mr =" << mr << endl;

if (m>mpr) {
ndv=ndv+1;
if (ndv>20) { //Lin Li: the value used to be 2,;
// sometimes not enough
cout <<"surface iteration did not converge" << endl;
stop;
}// if
}
else{
ndv=0;
}// if

//2011-05-18 Replaced goto 100 statement
if(m<=0) exit D100;
}// do D100;

if (n2>ibmx) {
cout <<"ibnd upper bound " << n2 << " exceeds ibmx" << endl;
stop;
}//if

//2011-05-18 Memory cleanup is done with DEALLOCATE F95
//statement
if(allocated(cbn1)) deallocate(cbn1);
if(allocated(cbn2)) deallocate(cbn2);
if(allocated(cbal)) deallocate(cbal);

if (bVerbose) {
cout <<"no. cavity mid-points inaccessible to solvent = " << ncav << endl;
}// if

//consolidate the list, removing dead boundary points, adding
//new ones..
j=0;
ncms=0;

//2011-05-19 Array is re-sized keeping old values in the
//memory.
if (allocated(ibgrd)) {
Nar=size(ibgrd);
if (Nar<ibmx) {
allocate(ibgrd_temp(Nar));
ibgrd_temp=ibgrd;
deallocate(ibgrd);
allocate(ibgrd(ibmx));
ibgrd(1:Nar)=ibgrd_temp;
deallocate(ibgrd_temp);
}// if
}
else{
allocate(ibgrd(ibmx));
}// if

for(i=1;i<=n2;i++){
ixyz=ibnd(i);
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;

if (bndeps[ix][iy][iz][1]!=0) {
j=j+1;
bndeps[ix][iy][iz][1]=j;

//2011-05-19 Precaution not to exceed array size (see
//above comment)
if (j<=ibmx) {
ibgrd(j)=ixyz;
}
else{
cout << "j=" << j << " is larger than ibmx= " << ibmx << " << thus stopped..." << endl;
stop;
}// if
}// if

for(jj=1;jj<=6;jj++){
//2011-05-17 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
ntt=comp(iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]]][imap[4][jj]);
nt=mod(ntt,epsdim);

if (nt<0) {
ntt=-ntt;

switch (imap[4][jj]){
case 1;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nX=ntt;
case 2;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nY=ntt;
case 3;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nZ=ntt;
case default;
cout <<"?????" << endl;
}// select;

if (nt==-1) {
ntt=ntt-1;

switch (imap[4][jj]){
case 1;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nX=ntt;
case 2;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nY=ntt;
case 3;
iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]].nZ=ntt;
case default;
cout <<"?????" << endl;
}// select;
}// if
}// if

if (bDebug) {
ntt=comp(iEpsMap[ix+imap[1][jj]][iy+imap[2][jj]][iz+imap[3][jj]]][imap[4][jj]);
nt=mod(ntt,epsdim);
jjj=ntt/epsdim;
if (nt==0&&jjj>0) {
cout <<"PROBLEMS 5" << ix << iy << iz << jj << endl;
}// if
}// if
}// do
}// do

if (j>ibmx){
cout << "no. ms points exceeds ibmx" << endl;
stop;
}// if

iBoundNum=j;

if (bVerbose) {
cout <<"after surface elaboration iBoundNum= " << iBoundNum << endl;
cout <<" and iBoundNumsurf= " << iBoundNumsurf << endl;
}// if

}// if

if (bDebug) {
open(52,file='bgpDurante.log',form='formatted');
write(52,*)iGrid,iBoundNum;
for(iiii=1;iiii<=iBoundNum;iiii++){
write(52,*)iiii,ibnd(iiii);
write(52,*)iiii,ibnd(3*iiii-2),ibnd(3*iiii-1),ibnd(3*iiii);
}// do
close (52);
}//if

if(allocated(bndeps)) deallocate(bndeps);
if(allocated(ibnd)) deallocate(ibnd);

//fScale bondary grid point positions relative to acc data
if (isolv&&(irea||logs||lognl||isen||isch)) {
if(bVerbose)cout <<"scaling boundary grid points" << endl;

//2011-05-19 Arrays allocated by ordinary F95 allocate
//statement
allocate(scspos(iBoundNum));

for(j=1;j<=iBoundNum;j++){
//2011-05-19 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
scspos(j)=Int2Float(ibgrd(j));
}// do

allocate(scsnor(iBoundNum),atsurf(iBoundNum),atndx(iBoundNum));

//2011-05-19 Other parameters to void are transfered via
//module architecture (declared in pointers and qlog modules)
call sclbp(scspos,scsnor,iBoundNum,iab1,iab2);

if (bVerbose) {
cout << iall << " points had to be assigned by global comparison" << endl;
}// if

if (!isite&&allocated(scsnor)) deallocate(scsnor);
}// if

if (isrf&&!ivertnorm){
if (bOnlyMol) {
//2011-05-19 Parameters transfered via module architecture
allocate(egrid(iGrid,iGrid,iGrid));

call msrf;

if(allocated(egrid) )deallocate(egrid);
}
else{
cout <<"msrf routine cannot be run" << endl;
cout <<"because there are also geometric objects" << endl;
}// if
}// if

if (!isitsf&&!isite&&!(isch&&scrgfrm!=0)) {
if(allocated(atndx)) deallocate(atndx);
if(allocated(atsurf)) deallocate(atsurf);
}// if

if(allocated(iab1)) deallocate(iab1);
if(allocated(iab2)) deallocate(iab2);
if(allocated(icume)) deallocate(icume);
if(allocated(r0)) deallocate(r0);
if(allocated(r02)) deallocate(r02);
if(allocated(rs2)) deallocate(rs2);
if(allocated(ast)) deallocate(ast);

}// void vwtms;
