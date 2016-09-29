#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;
// ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//==================================================================
// 2011-05-12 Parameters transfered via modules
void setiGaussian(radprobe);

// void setout(iNatom,radprobe,fExternRadius,fScale,iGrid,iNObject
// & ,cOldMid,iNMedia,bOnlyMol,isitmd,bDebug,bVerbose)
// include 'pointer.h'

// dimension xn2[3][iNatom]][rad3(iNatom)][xn(3)
// dimension iepsmp(iGrid,iGrid,iGrid,3)
//!c b+++++++++++
// int iNObject,imedia,objecttype,epsdim,kind,iac,iNMedia
// character*96 dataobject[iNObject][2]][strtmp][strtmp1
// dimension iAtomMed[iNatom+iNObject]
// float sLimGridUnit[3][2][iNObject]][vectz(3)
// float radprobe,exrd,modul2,axdist
// int ioff[3][1]][ismin(3)][ismax(3)
// dimension iepsmp(iGrid,iGrid,iGrid,3)
// int iNObject,imedia,objecttype,epsdim,kind,iac,iNMedia
// character*96 dataobject[iNObject][2]][strtmp][strtmp1
// dimension iAtomMed[iNatom+iNObject]
// float sLimGridUnit[3][2][iNObject]][vectz(3)
// float radprobe,exrd,modul2,axdist
// float fRMid,cOldMid(3),tmpvect(3),tmp,tmp1,dx,dy,dz,dist,zeta
// float xa(3),xb(3),radius,modul,x2,y2,z2,tmpvect1(3),mod2,tmp2
// float tmpvect2(3),xc(3),xd(3),xp(3),alpha,tan2,dot,modx,mody
// float vectx(3),vecty(3),dx1,dx2,dx3,dis2min1,dis2min2
// logical*1 bDebMap[iGrid][iGrid][iGrid]
// bool isitmd,bVerbose ! qlog module
// 2011-05-12 Some local variables
// float sq1(-15:15),sq2(-15:15),sq3(-15:15)
// float rad2a1(-15:15),rad2a2(-15:15),rad2a3(-15:15)
// SGrid <float> sq[-15:15] ][ rad2aav[-15:15]
SGrid <float> sq[-150:150],rad2aav[-150:150],vtemp2[-150:150]//Gaussian enlarged the values;
bool itobig,itest2,ipore,bOnlyMolbDebug;
character(96) :: strtmp,strtmp1;
// 2011-05-12 Non-standard int variable, thus necessary
int epsdim, objecttype,iflag;
// 2011-05-12 Non-standard float variable, thus necessary
float modul,modul2, mod2,modx,mody,dentemp;
// 2011-05-12 Non-standard type of variable, thus necessary
SGrid <float> xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist;
SGrid <float> tmpvect1,tmpvect2,origin;
SGrid <float> vectx,vecty,vectz,rad2av,fxn,vtemp;
SGrid <int> ismin,ismax,idist,idist1,ixyz,itest,ixn,i123;
real*8 ::coeff,stepsize;
integer*8 ::longint;
//!c b+++++++++++
//!c here fRadPrb is not zero only if one wants to map the ext}//ed surf.
//!c that has to be done with fRadPrb[1] if solute-solvent interface is
//!c concerned
//!c imedia = medium number for a object

//!c a non-zero entry in iepsmp indicates an atom # plus 1 (to properly treat
//!c re-entrant mid-points later (15 Aug 93)

//------------------------------------------------------------------------------------------
// 2011-05-27 Declarations added due to IMPLICIT NONE
int iboxt,iac,ibox,ii,igrdc,i,j,k,imedia,iv,ix,iy,iz,kind;
int limmax,lim,n;
integer,dimension(1:6) ::inwater //mid point;
logical,dimension(1:6) ::ifinwater //mid point;
float alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2;
float rad2a,rad4,radius,radmax2,fRMid,rad5,radp2,radtest,radprobe;
float radp,tan2,temp,tmp,tmp1,tmp2,epstemp;


//------------------------------------------------------------------------------------------

if(iGaussian==1&&inhomo==0&&logs) {
goto 1010;
}//if

pi=3.14159265359;
sigmatime=3.0;

gepsmp2=coord(repsout,repsout,repsout);
gepsmp=coord(0.0,0.0,0.0);
if(bVerbose)cout << "Starting creating Van der Waals Epsilon Map " << endl;
//!c e+++++++++++

epsdim=iNatom+iNObject+2;
//bDebug cout <<"-----setout-->" << epsdim << size(iAtomMed)
iboxt=0;
radmax2=0.0;
fRMid=Int2Float((iGrid+1)/2);
itest2=false;
if (iNatom>0) {
// 2011-05-12 Label 691 removed
for(ix=1;ix<=iNatom;ix++){
// 2011-05-13 Chaged to derive-type array sDelPhiPDB (pointers module)
radmax2=max(radmax2][sDelPhiPDB[ix].radius);
}// do
//691 continue
//!c b+++++++++++++++++++++++++++++
//!c this is probably the best way to do it,dep}//ing on which surf. is desired
temp=max(radprobe,fExternRadius);
// cout << << "Gaussian:in setout:temp << radprobe << fExternRadius << radmax2" << temp << radprobe << fExternRadius << radmax2
//!c e++++++++++++++++++++++++++++++
// radmax2=fScale*(radmax2+temp)
radmax2=sigmatime*fScale*(radmax2*sigma+temp) //Gaussian: 3 sigma plus temp. sigma=sigma*radius.Now radmax is 3 sigma.;
lim=1+radmax2;
// cout << << "Gaussian:setout:radmax2 << lim" << radmax2 << lim
limmax = 1200 //Gaussian:original value:12;
itobig=false;
if(lim>limmax) itobig=true;
igrdc=(2*lim+1)**3;
allocate(ioff[igrdc]);
if(!itobig) {
// 2011-05-12 removed goto 7878 statement
// if(itobig) goto 7878
radtest= (radmax2 + 0.5*sqrt(3.0))**2;
ibox=0;
// print* << "Gaussian:temp << radmax2 << " << temp << radmax2
// 2011-05-12 Strange statement. May allocate or may not allocate array that
// used later in the program irrespectively of itobig value, thus
// moved array allocation before if condition
// i_ioff= memalloc(i_ioff,4,3*igrdc)
for(ix=-lim;ix<=lim;ix++){
for(iy=-lim;iy<=lim;iy++){
for(iz=-lim;iz<=lim;iz++){
idist=int_coord(ix,iy,iz);
dist=Int2Float(idist.dot.idist);
ddist = dist + 0.25 + Int2Float(idist);
if((dist<radtest)||(ddist.vorlt.radtest)) {
ibox=ibox+1;
ioff[ibox]=idist;
}
else{
}// if
}// do
}// do
}// do
}// if
// 2011-05-12 Labels replaced by }// do and end if.
//! c +++++++++++++++++++++++++
}// if
//!c set interiors in OBJECTS
for(ii=1;ii<=iNObject;ii++){
ipore=false;
strtmp=dataobject[ii][1];
cout <<" " << endl;
if (strtmp(1:4)!='is a') {
strtmp1=dataobject[ii][2];
read(strtmp(8:10),'(I3)')objecttype;
read(strtmp(12:14),'(I3)')imedia;

//!c completing iAtomMed with imedia data
iAtomMed[iNatom+ii]=imedia;
if(bVerbose)cout << "(setout) object number" << ii << " in medium" << imedia << endl;
//!c check if object sits within limits of box and calculating ismin
//!c and ismax accordingly
exrd=fExternRadius*fScale;
temp=radprobe*fScale+exrd;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ismin= Float2Int(sLimGridUnit[ii].nMin-temp-1.);
ismin=min(ismin,iGrid);
ismin=max(ismin,1);
ismax= Float2Int(sLimGridUnit[ii].nMax+temp+1.);
ismax=min(ismax,iGrid);
ismax=max(ismax,1);

// 2011-05-13 Changed multiple IF's to switch
switch (objecttype){
case 1// if (objecttype==1) {
//!c dealing with a sphere
read(strtmp(16:80),*)kind,xb,radius;
ipore=(kind==2);
radius=radius*fScale;
radp=radius+exrd;
rad2=radius*radius;
radp2=radp*radp;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
xb=(xb-cOldMid)*fScale+fRMid;
for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
ixyz=int_coord(ix,iy,iz);
xyz2=(Int2Float(ixyz)-xb)*(float(ixyz)-xb);
if(sum(xyz2)<=radp2) bDebMap[ix][iy][iz]=ipore;
xyz2.nX=xyz2.nX+0.25;

ddist=sum(xyz2)+Int2Float(ixyz)-xb;

if(ddist.nX<=rad2) iepsmp(ix,iy,iz).nX=iNatom+1+ii+imedia*epsdim;
if(ddist.nY<=rad2) iepsmp(ix,iy,iz).nY=iNatom+1+ii+imedia*epsdim;
if(ddist.nZ<=rad2) iepsmp(ix,iy,iz).nZ=iNatom+1+ii+imedia*epsdim;
}// do
}// do
}// do


case 2 //if (objecttype==2) {
//!c dealing with a cylinder
read(strtmp(16:80),*)kind,xa,xb,radius;
ipore=(kind==2);
read(strtmp1,'(5f8.3)')vectz,modul,modul2;
//!c here we work in grid units
radius=radius*fScale;
rad2=radius*radius;
radp=radius+exrd;
radp2=radp*radp;
modul=modul*fScale;
modul2=modul*modul;
tmp=exrd*modul;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
xb=(xb-cOldMid)*fScale + fRMid;
vectz=vectz*fScale;
for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
//!c vectz=A-B
modul=|A-B|;
tmpvect1=P-B;
tmp1=(A-B)(P-B);
//!c mod2=(P-B)**2
modul2=(A-B)**2, tmp=fExternRadius*|A-B|.;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ixyz=int_coord(ix,iy,iz);
tmpvect1=Int2Float(ixyz)-xb;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=-tmp)&&(tmp1<=modul2+tmp)) {
mod2=tmpvect1.dot.tmpvect1;
if((mod2-(tmp1/modul)**2)<=radp2) bDebMap[ix][iy][iz]=ipore;
}// if

// 2011-05-13 x-offset
tmpvect1.nX=tmpvect1.nX+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iepsmp(ix,iy,iz).nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if
// 2011-05-13 y-offset
tmpvect1.nX=tmpvect1.nX-.5;
tmpvect1.nY=tmpvect1.nY+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iepsmp(ix,iy,iz).nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if
// 2011-05-13 z-offset
tmpvect1.nY=tmpvect1.nY-.5;
tmpvect1.nZ=tmpvect1.nZ+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iepsmp(ix,iy,iz).nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// do
}// do
}// do

case 3 // if (objecttype==3) {
//!c dealing with a cone
read(strtmp(16:80),*)kind,xa,xb,alpha;
ipore=(kind==2);
//!c conversion degrees --> radiants
alpha=alpha*3.14159265359/180.;
tan2=tan(alpha*.5)**2;
read(strtmp1,'(5f8.3)')vectz,modul,modul2;
//!c here we work in grid units
modul=modul*fScale;
modul2=modul*modul;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
xb=(xb-cOldMid)*fScale + fRMid;
vectz=vectz*fScale;
for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
ixyz=int_coord(ix,iy,iz);
tmpvect1=(Int2Float(ixyz)-fRMid)/fScale+cOldMid;
//!cWWW if ((ix==17)&&(iy==6)&&(iz==13)) {
//!c cout << ismin << ismax
//!c cout << tmpvect1
//!c cout << vectz << modul << modul2
//!c }// if
// 2011-05-13 Parameters transfered via module architecture
call distobj(tmpvect1,dist,vecdist,ii,fExternRadius,true);
if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;

//!c vectz=A-B
tmpvect1=P-B;
tmp1=(A-B)(P-B);
//!c mod2=(P-B)**2
modul2=(A-B)**2.;

// 2011-05-13 x-offset
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
tmpvect1=Int2Float(ixyz)-xb;
tmpvect1.nX=tmpvect1.nX+0.5;
tmp1=tmpvect1.dot.vectz;

if ((tmp1>=0.0)&&(tmp1<=modul2)) {
//!c mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)
//!c tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
mod2=(tmpvect1.dot.tmpvect1)*modul2;
// mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iepsmp(ix,iy,iz).nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if
// 2011-05-13 y-offset
tmpvect1.nX=tmpvect1.nX-0.5;
tmpvect1.nY=tmpvect1.nY+0.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
//!c mod2=tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2
//!c tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
mod2=(tmpvect1.dot.tmpvect1)*modul2;
// mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iepsmp(ix,iy,iz).nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if
// 2011-05-13 z-offset
tmpvect1.nY=tmpvect1.nY-0.5;
tmpvect1.nZ=tmpvect1.nZ+0.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
//!c mod2=tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2
//!c tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
mod2=(tmpvect1.dot.tmpvect1)*modul2;
// mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iepsmp(ix,iy,iz).nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if

}// do
}// do
}// do
// }// if

case 4 // if (objecttype==4) {
//!c dealing with a parallelepiped
read(strtmp(16:80),*)kind,xa,xb,xc,xd;
ipore=(kind==2);
read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx, modx;

//!c conversion to axial symmetry points
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
xb=0.5*(xc+xb);
xb=(xb-cOldMid)*fScale + fRMid;

modul=modul*fScale;
modul2=modul*modul;
vectz=vectz*fScale;
modx=modx*fScale;
tmp1=modx*modx/2.;
vectx=vectx*fScale;
mody=mody*fScale;
tmp2=mody*mody/2.;
vecty=vecty*fScale;
for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
//!c vectz=A-B
vectx=C-D;
vecty=(C+D)/2-B;
tmp1=|C-D|/2;
//!c modul2=(A-B)**2, tmp2=mody/2
//!c dot=(P-B)(..-..)

ixyz=int_coord(ix,iy,iz);
xp=(Int2Float(ixyz)-fRMid)/fScale + cOldMid;
// 2011-05-13 Parameters transfered via module architecture
call distobj(xp,dist,vecdist,ii,fExternRadius,true);
if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;
//!c now xp=P-B

// 2011-05-13 x-offset
xp=Int2Float(ixyz)-xb;
xp.nX=xp.nX+0.5;
dot=vectz.dot.xp;
if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
if (abs(dot)<=tmp2) {
iepsmp(ix,iy,iz).nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if
// 2011-05-13 y-offset
xp.nX=xp.nX-0.5;
xp.nY=xp.nY+0.5;
dot=vectz.dot.xp;
if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
if (abs(dot)<=tmp2) {
iepsmp(ix,iy,iz).nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if
// 2011-05-13 z-offset
xp.nY=xp.nY-0.5;
xp.nZ=xp.nZ+0.5;
dot=vectz.dot.xp;
// xp(1)=ix-xb(1)
// xp(2)=iy-xb(2)
// xp(3)=iz+.5-xb(3)
// dot=vectz(1)*xp(1)+vectz(2)*xp(2)+vectz(3)*xp(3)
if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
// dot=vectx(1)*xp(1)+vectx(2)*xp(2)+vectx(3)*xp(3)
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
// dot=vecty(1)*xp(1)+vecty(2)*xp(2)+vecty(3)*xp(3)
if (abs(dot)<=tmp2) {
iepsmp(ix,iy,iz).nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if

}// do
}// do
}// do
}// select;

}// if
}// do

//!c }// of setting in OBJECTS
//!c e++++++++++++++++++++++++

//!c set interiors in MOLECULES

if(itest2||itobig) cout <<"setout method 1" << itest2 << itobig << endl;
DoATOMS: do iv=1, iNatom;
//!c restore values
rad= sDelPhiPDB[iv].radius;
rad= sigmatime*sigma*rad // 3 sigma;

// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
xn=xn2[iv];
// xn(1)=xn2[1][iv]
xn(2)=xn2[2][iv];
xn(3)=xn2[3][iv];
// 2011-05-13 Removed GOTO statement

if(rad<1.e-6) {
//bDebug cout <<"------setout22:iv-->" << iv << ixyz << iepsmp(ix << iy << iz) << epsdim << iAtomMed[iv]
cycle DoATOMS;
}// if
// if(rad<1.e-6) goto 608
//!c fScale radius to grid
// cout << << "rad << fScale << rad*fScale << radprobe" << rad << fScale << rad*fScale << radprobe
rad= rad*fScale;
rad5= (rad + 0.5)**2;
radp = rad + fExternRadius*fScale;
//!c b++++++++++++++++++++++++++++
rad= rad + radprobe*fScale;
//!c e++++++++++++++++++++++++++++rad2 is sigmatimed rad square. radsq is original rad square.
rad4= (rad + 0.5)**2;
rad2= rad*rad;
radp2= radp*radp;
radsq=rad2/sigmatime**2;


//!c set dielectric map
//!c check if sphere sits within limits of box
itest2=false;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ismin=Float2Int(xn-radmax2-1.);
ismax=Float2Int(xn+radmax2+1.);
itest=ismin;
ismin=min(ismin,iGrid);
ismin=max(ismin,1);
if(itest.vorne.ismin) itest2=true;
itest=ismax;
ismax=min(ismax,iGrid);
ismax=max(ismax,1);
if(itest.vorne.ismax) itest2=true;
//--------------Maxim Sep_15, 2011---------------------------
// if(itest2) cout <<"-->" << iv << xn << sDelPhiPDB[iv].radius << ismin << radmax2
//-----------------------------------------------------------


//!--------- slow method--------------------------
// if(itest2||itobig) {
if(itobig) {
//bDebug cout <<"==itest2====>" << iv << iNMedia << bOnlyMol

// cout << << "setout:slow method"
// cout << << "itest2 << itobig" << itest2 << itobig
// 2011-05-13 Seems redundant statement
// num=num+1
rad2a = rad2 - 0.25;
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(ix=ismin.nX;ix<=ismax.nX;ix++){
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ixyz=int_coord(ix,iy,iz);
dxyz=Int2Float(ixyz)-xn;
distsq=dxyz.dot.dxyz;
dxyz=dxyz+distsq;

//!c b+++++++++++++++++++
if(dxyz.nX<rad2a) {
iepsmp(ix][iy][iz).nX=iv+1+iAtomMed[iv]*epsdim;
}// if
if(dxyz.nY<rad2a) {
iepsmp(ix][iy][iz).nY=iv+1+iAtomMed[iv]*epsdim;
}// if
if(dxyz.nZ<rad2a) {
iepsmp(ix][iy][iz).nZ=iv+1+iAtomMed[iv]*epsdim;
}// if
//!c e+++++++++++++++++++
if(distsq<radp2) bDebMap[ix][iy][iz] =false;
}// do
}// do
}// do
// 2011-05-13 Labels removed
}
else{


//!c faster method
//!c IT HAS PROBLEMS!!!! Walter (be careful before using
//!c also in multidilectric case!!!&&!isitmd

rad2a = rad2 - 0.25;
// 2011-05-13 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ixn=nFloat2Int(xn);
fxn=Int2Float(ixn)-xn;
rad2av=rad2a-fxn //rad2av is never used;
for(ix=-lim;ix<=lim;ix++){
vtemp=Int2Float(ix)+fxn;
sq[ix]=vtemp*vtemp;
rad2aav[ix]=rad2a-vtemp;
vtemp2[ix]=vtemp;
}// do
//!c all_f2c.sh array_globle.lst array.lst auto_cpp backup f2c.sh f2c_v2.tmp f2c_v.tmp header.txt ijk_kji inc_epsmakmod_setout.cpp inc_epsmakmod_setout.f95 list.txt log.txt qdiff.cpp qdiff.f95 setout.cpp space_module_fortran space_module_fortran2 src test.sh test.txt variables.list vwtms.cpp
// cout << << "Gaussian:iNMedia << bOnlyMol:" << iNMedia << bOnlyMol
//!c b+++adjust inter-atom, different epsilon bgps+++04/2004 Walter+
if(iNMedia>1&&bOnlyMol) {
// cout << << "Gaussian:setout:fast 1 method"
for(i=1;i<=ibox;i++){
// do 2004 i=1,ibox
// 2011-05-14 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
i123=ioff[i];
ixyz=ixn+i123;
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;
distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
cout <<"------->" << i << ibox << i123 << ixyz << distsq << endl;

//bDebug cout <<"-1->" << iv << distsq << rad2aav[i123.nX].nX << rad2aav[i123.nY].nY << &
// & rad2aav[i123.nZ].nZ
if(distsq<rad2aav[i123.nX].nX) {
// if(distsq<rad2a1(i1)) {
iac=mod(iepsmp(ix,iy,iz).nX,epsdim)-1;
// iac=mod(iepsmp(ix,iy,iz,1),epsdim)-1
if(iac==-1||iac>iNatom) {
//!c occhio! non so cosa fa con i pori!!
iepsmp(ix][iy][iz).nX=iv+1+iAtomMed[iv]*epsdim;
}
else{
// 2011-05-14 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nX= ddxyz.nX+0.5;
dis2min1=ddxyz.dot.ddxyz - rad2;
// dx1=ix+0.5-xn(1)
// dx2=iy-xn(2)
// dx3=iz-xn(3)
// dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nX= ddxyz.nX+0.5;
dis2min2=ddxyz.dot.ddxyz -(sDelPhiPDB[iac].radius*fScale)**2;
// dx1=ix-xn2[1][iac]+0.5
// dx2=iy-xn2[2][iac]
// dx3=iz-xn2[3][iac]
// dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*fScale)**2

if (dis2min2>dis2min1) iac=iv;
iepsmp(ix][iy][iz).nX=iac+1+iAtomMed[iac]*epsdim;
}//if
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 1 << 2) << "asse 1"
}// if
if(distsq<rad2aav[i123.nY].nY) {
// if(distsq<rad2a2(i2)) {
iac=mod(iepsmp(ix,iy,iz).nY,epsdim)-1;
// iac=mod(iepsmp(ix,iy,iz,2),epsdim)-1
if(iac==-1||iac>iNatom) {
//!c occhio! non so cosa fa con i pori!!
iepsmp(ix][iy][iz).nY=iv+1+iAtomMed[iv]*epsdim;
// iepsmp(ix][iy][iz][2)=iv+1+iAtomMed[iv]*epsdim
}
else{
// 2011-05-14 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nY= ddxyz.nY+0.5;
dis2min1=ddxyz.dot.ddxyz - rad2;
// dx1=ix-xn(1)
// dx2=iy+0.5-xn(2)
// dx3=iz-xn(3)
// dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nY= ddxyz.nY+0.5;
dis2min2=ddxyz.dot.ddxyz -(sDelPhiPDB[iac].radius*fScale)**2;
// dx1=ix-xn2[1][iac]
// dx2=iy+0.5-xn2[2][iac]
// dx3=iz-xn2[3][iac]
// dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*fScale)**2
if (dis2min2>dis2min1) iac=iv;
iepsmp(ix][iy][iz).nY=iac+1+iAtomMed[iac]*epsdim;
// iepsmp(ix][iy][iz][2)=iac+1+iAtomMed[iac]*epsdim
}//if
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 2 << 2) << "asse 2"
}// if

if(distsq<rad2aav[i123.nZ].nZ) {
// if(distsq<rad2a3(i3)) {
iac=mod(iepsmp(ix,iy,iz).nZ,epsdim)-1;
// iac=mod(iepsmp(ix,iy,iz,3),epsdim)-1
if(iac==-1||iac>iNatom) {
//!c occhio! non so cosa fa con i pori!!
iepsmp(ix][iy][iz).nZ=iv+1+iAtomMed[iv]*epsdim;
// iepsmp(ix][iy][iz][3)=iv+1+iAtomMed[iv]*epsdim
}
else{
// 2011-05-14 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nZ= ddxyz.nZ+0.5;
dis2min1=ddxyz.dot.ddxyz - rad2;
// dx1=ix-xn(1)
// dx2=iy-xn(2)
// dx3=iz+0.5-xn(3)
// dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nZ= ddxyz.nZ+0.5;
dis2min2=ddxyz.dot.ddxyz -(sDelPhiPDB[iac].radius*fScale)**2;
// dx1=ix-xn2[1][iac]
// dx2=iy-xn2[2][iac]
// dx3=iz+0.5-xn2[3][iac]
// dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*fScale)**2
if (dis2min2>dis2min1) iac=iv;
iepsmp(ix][iy][iz).nZ=iac+1+iAtomMed[iac]*epsdim;
// iepsmp(ix][iy][iz][3)=iac+1+iAtomMed[iac]*epsdim
}//if
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 3 << 2) << "asse 3"
}// if

if(distsq<radp2) bDebMap[ix][iy][iz]=false;
}// do
//2004 continue
}
else{
//!c all_f2c.sh array_globle.lst array.lst auto_cpp backup f2c.sh f2c_v2.tmp f2c_v.tmp header.txt ijk_kji inc_epsmakmod_setout.cpp inc_epsmakmod_setout.f95 list.txt log.txt qdiff.cpp qdiff.f95 setout.cpp space_module_fortran space_module_fortran2 src test.sh test.txt variables.list vwtms.cpp




// cout << << "Gaussian:setout:fast 2 method"
//!C$DIR NO_RECURRENCE
// cout <<"Gaussian:iv << iNMedia << bOnlyMol << ibox:" << iv << iNMedia << bOnlyMol << ibox
for(i=1;i<=ibox;i++){


// do 9024 i=1,ibox
// 2011-05-14 Using operations on coord and int_coord type variables defined
// in module operators_on_coordinates
i123=ioff[i];
ixyz=ixn+i123;
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;
//Gaussian:Avoid reaching grids out side (1,iGrid)
if(ix<iGrid&&ix>1&&iy<iGrid&& iy>1&&iz<iGrid&&iz>1) {
distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
// i1= ioff[1][i]
i2= ioff[2][i];
i3= ioff[3][i];
// ix=ixn1+ i1
iy=iyn1+ i2;
iz=izn1+ i3;
// distsq = sq1(i1) +sq2(i2) + sq3(i3)
//!c b+++++++++++++++++++
//bDebug cout <<"-2->" << iv << distsq << rad2aav[i123.nX].nX << rad2aav[i123.nY].nY << &
// & rad2aav[i123.nZ].nZ
// cout <<"-2a->" << i << ibox << i123 << ixn
// cout <<"-2b->" << i << ibox << sq[i123.nX].nX << sq[i123.nY].nY << sq[i123.nZ].nZ
// cout << << "Gaussian:distsq << rad2aav[i123.nX].nX:" << distsq << rad2aav[i123.nX].nX

if(distsq<rad2aav[i123.nX].nX) {
// iepsmp(ix][iy][iz).nX=iv+1+iAtomMed[iv]*epsdim
// gepsmp(ix,iy,iz).nX=repsin
distance2=distsq+0.25+vtemp2[i123.nX].nX;
// den=peak*exp(-(distance2/(sigma**2*radsq)))
den=exp(-(distance2/(sigma**2*radsq)));
gepsmp(ix,iy,iz).nX=1-(1-gepsmp(ix,iy,iz).nX)*(1-den);

// gepsmp2(ix,iy,iz).nX=gepsmp(ix,iy,iz).nX*repsin+(1-gepsmp(ix,iy,iz).nX)*repsout
// if(distsq<rad2a1(i1)) {
// iepsmp(ix][iy][iz][1)=iv+1+iAtomMed[iv]*epsdim
// cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz).nX << "asse 1"
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 1 << 2) << "asse 1"
}// if

if(distsq<rad2aav[i123.nY].nY) {
// iepsmp(ix][iy][iz).nY=iv+1+iAtomMed[iv]*epsdim
// gepsmp(ix,iy,iz).nY=repsin
distance2=distsq+0.25+vtemp2[i123.nY].nY;
den=exp(-(distance2/(sigma**2*radsq)));

// if(ix==100&&iy==100&&iz==100) &
// den=peak*exp(-(distance2/(sigma**2*radsq)))
gepsmp(ix,iy,iz).nY=1-(1-gepsmp(ix,iy,iz).nY)*(1-den);
// gepsmp2(ix,iy,iz).nY=gepsmp(ix,iy,iz).nY*repsin+(1-gepsmp(ix,iy,iz).nY)*repsout

// iepsmp(ix][iy][iz][2)=iv+1+iAtomMed[iv]*epsdim
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 2 << 2) << "asse 2"
}// if

if(distsq<rad2aav[i123.nZ].nZ) {
// iepsmp(ix][iy][iz).nZ=iv+1+iAtomMed[iv]*epsdim
// gepsmp(ix,iy,iz).nZ=repsin
distance2=distsq+0.25+vtemp2[i123.nZ].nZ;
den=exp(-(distance2/(sigma**2*radsq))) //Gaussian: make density at center is 1.;
// den=peak*exp(-(distance2/(sigma**2*radsq)))
gepsmp(ix,iy,iz).nZ=1-(1-gepsmp(ix,iy,iz).nZ)*(1-den);
// gepsmp2(ix,iy,iz).nZ=gepsmp(ix,iy,iz).nZ*repsin+(1-gepsmp(ix,iy,iz).nZ)*repsout
// if(distsq<rad2a1(i1)) {
//
// if(distsq<rad2a3(i3)) {
// iepsmp(ix][iy][iz][3)=iv+1+iAtomMed[iv]*epsdim
//!c cout << "atomo:" << iv << " mezzo:" << iepsmp(ix << iy << iz << 3 << 2) << "asse 3"
}// if

//!c e++++++++++++++++++++
// if(distsq<radp2) bDebMap[ix][iy][iz]=false
//Gaussian change method of generating bDebMap

}//if //for detecting outside cube;
}// do
//!9024 continue
}// if
}// if
// 2011-05-13 Removed label 608
// 608 continue
// cout <<"------setout22:iv-->" << iv << ixyz << iepsmp(ix << iy << iz) << epsdim << iAtomMed[iv]

//!c }// do of atoms

}// do DoATOMS;


// cutoff=0.9
// cout << << "cutoff << sigma:" << cutoff << sigma
if(false){ //cutoff: make the top flat;
for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
for(iz=1;iz<=iGrid;iz++){
gepsmp(ix,iy,iz).nX=gepsmp(ix,iy,iz).nX/cutoff;
gepsmp(ix,iy,iz).nY=gepsmp(ix,iy,iz).nY/cutoff;
gepsmp(ix,iy,iz).nZ=gepsmp(ix,iy,iz).nZ/cutoff;
if(gepsmp(ix,iy,iz).nX>1.0) gepsmp(ix,iy,iz).nX=1.0;
if(gepsmp(ix,iy,iz).nY>1.0) gepsmp(ix,iy,iz).nY=1.0;
if(gepsmp(ix,iy,iz).nZ>1.0) gepsmp(ix,iy,iz).nZ=1.0;
}//do
}//do
}//do
}//if

1010 continue;
// cout << << "bDebug: srfcut << inhomo:" << srfcut << inhomo
for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
for(iz=1;iz<=iGrid;iz++){
gepsmp2(ix,iy,iz).nX=gepsmp(ix,iy,iz).nX*repsin+(1-gepsmp(ix,iy,iz).nX)*repsout;
gepsmp2(ix,iy,iz).nY=gepsmp(ix,iy,iz).nY*repsin+(1-gepsmp(ix,iy,iz).nY)*repsout;
gepsmp2(ix,iy,iz).nZ=gepsmp(ix,iy,iz).nZ*repsin+(1-gepsmp(ix,iy,iz).nZ)*repsout;
//################### for set epsout in protein larger than in water:##########
if(gepsmp(ix,iy,iz).nX<0.02){
gepsmp2(ix,iy,iz).nX=80.0;
}//if
if(gepsmp(ix,iy,iz).nY<0.02){
gepsmp2(ix,iy,iz).nY=80.0;
}//if
if(gepsmp(ix,iy,iz).nZ<0.02){
gepsmp2(ix,iy,iz).nZ=80.0;
}//if
//################### }// for this epsout in protein different than in water######
}//do
}//do
}//do


if(inhomo==1){ //reduce epsilon out side protein;
// epstemp=repsin*(1-srfcut)+repsout*srfcut
epstemp=srfcut;
if(epstemp<repsin) {
print* << "srfcut is lower than epsin." << endl;
stop;
}//if
// print* << "epstemp" << epstemp
for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
for(iz=1;iz<=iGrid;iz++){
if(gepsmp2(ix,iy,iz).nX>epstemp){
gepsmp2(ix,iy,iz).nX=repsin;
}// if
if(gepsmp2(ix,iy,iz).nY>epstemp){
gepsmp2(ix,iy,iz).nY=repsin;
}// if
if(gepsmp2(ix,iy,iz).nZ>epstemp){
gepsmp2(ix,iy,iz).nZ=repsin;
}// if
}//do
}//do
}//do
}//if

iBoundNum=0;
longint=0;
dentemp=0.1;
// ifinwater=true
for(i=2;i<=iGrid;i++){
for(j=2;j<=iGrid;j++){
for(k=2;k<=iGrid;k++){
if(gepsmp(i,j,k).nX>dentemp||gepsmp(i,j,k).nY>dentemp & ||gepsmp(i,j,k).nZ>dentemp||gepsmp(i-1,j,k).nX>dentemp ||gepsmp(i,j-1,k).nY>dentemp||gepsmp(i,j,k-1).nZ>dentemp){

longint=longint+1;
bDebMap[i][j][k]=false //iGaussian change method of generating bDebMap;
}//if
}//do
}//do
}//do
iBoundNum=longint;
allocate(ibgrd[iBoundNum]);

n=0;
for(i=2;i<=iGrid;i++){
for(j=2;j<=iGrid;j++){
for(k=2;k<=iGrid;k++){

if(gepsmp(i,j,k).nX>dentemp||gepsmp(i,j,k).nY>dentemp & ||gepsmp(i,j,k).nZ>dentemp||gepsmp(i-1,j,k).nX>dentemp ||gepsmp(i,j-1,k).nY>dentemp||gepsmp(i,j,k-1).nZ>dentemp){


n=n+1;
ibgrd[n].nX=i;
ibgrd[n].nY=j;
ibgrd[n].nZ=k;
}//if
}//do
}//do
}//do

// Array deallocation is made by ordinary F95 statement
if(allocated(ioff)) deallocate(ioff);
// i_ioff= memalloc(i_ioff,0,0)
//!c b+++++bDebug+++++++++++++++++++++++++++++++
if (false) { //Gaussian;
open(52,file='linlieps_k',form='formatted');
write (52,*)iGrid,iGrid*iGrid;
for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
for(iz=1;iz<=iGrid;iz++){
write (52,*)ix,iy,iz,iepsmp(ix,iy,iz).nZ;
}// do
}// do
}// do
close (52);

}// if


if (false) { //bDebug: output density and eps;
open(14,file="cube_density", form='formatted');

write(14,*)'qdiffxs4 with an improved surfacing routine';
write(14,*) 'Gaussian cube format phimap';
coeff=0.5291772108;
stepsize=1.0/fScale;
origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;
// do i=1,3
// origin(i)=(cOldMid(i)-stepsize*(iGrid-1)/2)/coeff
// }//do
write(14,'(i5,3f12.6)') 1, origin;
// write(14,'(i5,3f12.6)') 1, (origin(i),i=1,3)
write(14,'(i5,3f12.6)') iGrid, stepsize/coeff,0.0,0.0;
write(14,'(i5,3f12.6)') iGrid, 0.0,stepsize/coeff,0.0;
write(14,'(i5,3f12.6)') iGrid, 0.0,0.0,stepsize/coeff;
write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0;

for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
// do iz=1,iGrid
write(14,'(6E13.5)')(gepsmp(ix,iy,iz).nZ,iz=1,iGrid);
// write(14,*),ix,iy,iz,gepsmp(ix,iy,iz).nX,gepsmp(ix,iy,iz).nY,gepsmp(ix,iy,iz).nZ

}// do
}// do
close(14);

open(14,file="cube_eps", form='formatted');

write(14,*)'qdiffxs4 with an improved surfacing routine';
write(14,*) 'Gaussian cube format phimap';
coeff=0.5291772108;
stepsize=1.0/fScale;
origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;
// do i=1,3
// origin(i)=(cOldMid(i)-stepsize*(iGrid-1)/2)/coeff
// }//do
write(14,'(i5,3f12.6)') 1, origin;
write(14,'(i5,3f12.6)') iGrid, stepsize/coeff,0.0,0.0;
write(14,'(i5,3f12.6)') iGrid, 0.0,stepsize/coeff,0.0;
write(14,'(i5,3f12.6)') iGrid, 0.0,0.0,stepsize/coeff;
write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0;

for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
write(14,'(6E13.5)')(gepsmp2(ix,iy,iz).nZ,iz=1,iGrid);
}// do
}// do
close(14);
}//if
//!c e++++++++++++++++++++++++++++++++++++++++++
if(bVerbose)cout <<"Ending creating Van der Waals Epsilon Map " << endl;

return;
}// void setiGaussian;

