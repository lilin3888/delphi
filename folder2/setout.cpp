#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>

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
// 2011-05-12 Parameters transfered via modules
//#####################################################################

void setout(radprobe);



//2011-05-12 Some local variables
SGrid <float> sq[-15:15] , rad2aav[-15:15];
bool itobig,itest2,ipore,bOnlyMolbDebug;
character(96) :: strtmp,strtmp1;
//2011-05-12 Non-standard int variable, thus necessary
int epsdim, objecttype;
//2011-05-12 Non-standard float variable, thus necessary
float modul,modul2, mod2,modx,mody;
//2011-05-12 Non-standard type of variable, thus necessary
SGrid <float> xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist;
SGrid <float> tmpvect1,tmpvect2;
SGrid <float> vectx,vecty,vectz,rad2av,fxn,vtemp;
SGrid <int> ismin,ismax,idist,idist1,ixyz,itest,ixn,i123;
//here fRadPrb is not zero only if one wants to map the ext}//ed
//surf. that has to be done with fRadPrb(1) if solute-solvent
//interface is concerned
//imedia = medium number for a object
//a non-zero entry in iEpsMap indicates an atom # plus 1 (to
//properly treat re-entrant mid-points later (15 Aug 93)
//2011-05-27 Declarations added due to IMPLICIT NONE
int iboxt,iac,ibox,ii,igrdc,i,imedia,iv,ix,iy,iz,kind;
int limmax,lim;
float alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2;
float rad2a,rad4,radius,radmax2,fRMid,rad5,radp2,radtest,radprobe;
float radp,tan2,temp,tmp,tmp1,tmp2;

if(bVerbose)cout << "Starting creating Van der Waals Epsilon Map " << endl;

epsdim=iNatom+iNObject+2;
iboxt=0;
radmax2=0.0;
fRMid=Int2Float((iGrid+1)/2);
itest2=false;

if (iNatom>0) {
//2011-05-12 Label 691 removed
for(ix=1;ix<=iNatom;ix++){
//2011-05-13 Chaged to derive-type array sDelPhiPDB (pointers
//module)
radmax2=max(radmax2][sDelPhiPDB[ix].radius);
}// do

//this is probably the best way to do it,dep}//ing on which
//surf. is desired
temp=max(radprobe,fExternRadius);

radmax2=fScale*(radmax2+temp);
lim=1+radmax2;
limmax = 12;
itobig=false;
if(lim>limmax) itobig=true;
igrdc=(2*lim+1)**3;
allocate(ioff[igrdc]);

if (!itobig) {
//2011-05-12 removed goto 7878 statement
radtest= (radmax2 + 0.5*sqrt(3.0))**2;
ibox=0;

//2011-05-12 Strange statement. May allocate or may not
//allocate array that used later in the program
//irrespectively of itobig value, thus moved array
//allocation before if condition
for(ix=-lim;ix<=lim;ix++){
for(iy=-lim;iy<=lim;iy++){
for(iz=-lim;iz<=lim;iz++){
//2011-05-13 Replaced by faster operation
//2011-05-13 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
idist=int_coord(ix,iy,iz);
dist=Int2Float(idist.dot.idist);
ddist = dist + 0.25 + Int2Float(idist);

if ((dist<radtest)||(ddist.vorlt.radtest)){
ibox=ibox+1;
ioff[ibox]=idist;
}// if
}// do
}// do
}// do
}// if
}// if

//set interiors in OBJECTS
for(ii=1;ii<=iNObject;ii++){
ipore=false;
strtmp=dataobject[ii][1];
cout <<" " << endl;

if (strtmp(1:4)!='is a') {
strtmp1=dataobject[ii][2];
read(strtmp(8:10),'(I3)')objecttype;
read(strtmp(12:14),'(I3)')imedia;

//completing iAtomMed with imedia data
iAtomMed[iNatom+ii]=imedia;
if(bVerbose)cout << "(setout) object number" << ii << " in medium" << imedia << endl;

//check if object sits within limits of box and calculating
//ismin and ismax accordingly
exrd=fExternRadius*fScale;
temp=radprobe*fScale+exrd;

//2011-05-13 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
ismin=Float2Int(sLimGridUnit[ii]%min-temp-1.);
ismin=min(ismin,iGrid);
ismin=max(ismin,1);

ismax=Float2Int(sLimGridUnit[ii]%max+temp+1.);
ismax=min(ismax,iGrid);
ismax=max(ismax,1);

//2011-05-13 Changed multiple IF's to switch
switch (objecttype){
case 1 //if (objecttype==1) {
//dealing with a sphere
read(strtmp(16:80),*)kind,xb,radius;
ipore=(kind==2);
radius=radius*fScale;
radp=radius+exrd;
rad2=radius*radius;
radp2=radp*radp;

//2011-05-13 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
xb=(xb-cOldMid)*fScale+fRMid;

for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
ixyz=int_coord(ix,iy,iz);
xyz2=(Int2Float(ixyz)-xb)*(float(ixyz)-xb);
if(sum(xyz2)<=radp2) bDebMap[ix][iy][iz]=ipore;
xyz2.nX=xyz2.nX+0.25;

ddist=sum(xyz2)+Int2Float(ixyz)-xb;

if(ddist.nX<=rad2) iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
if(ddist.nY<=rad2) iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
if(ddist.nZ<=rad2) iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
}// do
}// do
}// do
case 2 //if (objecttype==2) {
//dealing with a cylinder
read(strtmp(16:80),*)kind,xa,xb,radius;
ipore=(kind==2);
read(strtmp1,'(5f8.3)')vectz,modul,modul2;

//here we work in grid units
radius=radius*fScale;
rad2=radius*radius;
radp=radius+exrd;
radp2=radp*radp;
modul=modul*fScale;
modul2=modul*modul;
tmp=exrd*modul;

//2011-05-13 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
xb=(xb-cOldMid)*fScale + fRMid;
vectz=vectz*fScale;

for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
//vectz=A-B
modul=|A-B|;
tmpvect1=P-B;

//tmp1=(A-B)(P-B)
//mod2=(P-B)**2
modul2=(A-B)**2,;
//tmp=fExternRadius*|A-B|.
//2011-05-13 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
ixyz=int_coord(ix,iy,iz);
tmpvect1=Int2Float(ixyz)-xb;
tmp1=tmpvect1.dot.vectz;

if((tmp1>=-tmp)&&(tmp1<=modul2+tmp)){
mod2=tmpvect1.dot.tmpvect1;
if((mod2-(tmp1/modul)**2)<=radp2) bDebMap[ix][iy][iz]=ipore;
}// if

//x-offset
tmpvect1.nX=tmpvect1.nX+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if

//y-offset
tmpvect1.nX=tmpvect1.nX-.5;
tmpvect1.nY=tmpvect1.nY+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if

//z-offset
tmpvect1.nY=tmpvect1.nY-.5;
tmpvect1.nZ=tmpvect1.nZ+.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=tmpvect1.dot.tmpvect1;
if ((mod2-(tmp1/modul)**2)<=rad2) {
iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// do
}// do
}// do
case 3 //if (objecttype==3) {
//dealing with a cone
read(strtmp(16:80),*)kind,xa,xb,alpha;
ipore=(kind==2);

//conversion degrees --> radiants
alpha=alpha*3.14159265359/180.;
tan2=tan(alpha*.5)**2;
read(strtmp1,'(5f8.3)')vectz,modul,modul2;

//here we work in grid units
modul=modul*fScale;
modul2=modul*modul;

//2011-05-13 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
xb=(xb-cOldMid)*fScale + fRMid;
vectz=vectz*fScale;

for(ix=ismin.nX;ix<=ismax.nX;ix++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
ixyz=int_coord(ix,iy,iz);
tmpvect1=(Int2Float(ixyz)-fRMid)/fScale+cOldMid;

//2011-05-13 Parameters transfered via module
//architecture
call distobj(tmpvect1,dist,vecdist,ii,fExternRadius,true);

if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;

//vectz=A-B
tmpvect1=P-B;
tmp1=(A-B)(P-B);
//mod2=(P-B)**2
modul2=(A-B)**2.;

//x-offset
//2011-05-13 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
tmpvect1=Int2Float(ixyz)-xb;
tmpvect1.nX=tmpvect1.nX+0.5;
tmp1=tmpvect1.dot.vectz;

if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=(tmpvect1.dot.tmpvect1)*modul2;
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if

//y-offset
tmpvect1.nX=tmpvect1.nX-0.5;
tmpvect1.nY=tmpvect1.nY+0.5;
tmp1=tmpvect1.dot.vectz;

if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=(tmpvect1.dot.tmpvect1)*modul2;
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if

//z-offset
tmpvect1.nY=tmpvect1.nY-0.5;
tmpvect1.nZ=tmpvect1.nZ+0.5;
tmp1=tmpvect1.dot.vectz;
if ((tmp1>=0.0)&&(tmp1<=modul2)) {
mod2=(tmpvect1.dot.tmpvect1)*modul2;
tmp2=(1+tan2)*tmp1*tmp1;
if (mod2<=tmp2) {
iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if

}// do
}// do
}// do
case 4 //if (objecttype==4) {
//dealing with a parallelepiped
read(strtmp(16:80),*)kind,xa,xb,xc,xd;
ipore=(kind==2);
read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx,modx;

//conversion to axial symmetry points
//2011-05-13 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
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
//vectz=A-B
vectx=C-D;
vecty=(C+D)/2-B;

//tmp1=|C-D|/2
//modul2=(A-B)**2, tmp2=mody/2
//dot=(P-B)(..-..)

ixyz=int_coord(ix,iy,iz);
xp=(Int2Float(ixyz)-fRMid)/fScale + cOldMid;

//2011-05-13 Parameters transfered via module
//architecture
call distobj(xp,dist,vecdist,ii,fExternRadius,true);

if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;

//now xp=P-B

//x-offset
xp=Int2Float(ixyz)-xb;
xp.nX=xp.nX+0.5;
dot=vectz.dot.xp;

if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
if (abs(dot)<=tmp2) {
iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if

//y-offset
xp.nX=xp.nX-0.5;
xp.nY=xp.nY+0.5;
dot=vectz.dot.xp;

if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
if (abs(dot)<=tmp2) {
iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if

//z-offset
xp.nY=xp.nY-0.5;
xp.nZ=xp.nZ+0.5;
dot=vectz.dot.xp;

if ((dot>=0.0)&&(dot<=modul2)) {
dot=vectx.dot.xp;
if (abs(dot)<=tmp1) {
dot=vecty.dot.xp;
if (abs(dot)<=tmp2) {
iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
}// if
}// if
}// if

}// do
}// do
}// do
}// select;
}// if
}// do //end of setting in OBJECTS;

//set interiors in MOLECULES
if(itest2||itobig) cout <<"setout method 1" << itest2 << itobig << endl;

DoATOMS: do iv=1,iNatom;
//restore values
rad= sDelPhiPDB[iv].radius;

//Using operations on coord and int_coord type variables defined
//in module operators_on_coordinates
xn=xn2[iv];

//2011-05-13 Removed GOTO statement
if (rad<1.e-6) {
cycle DoATOMS;
}// if

//fScale radius to grid
rad=rad*fScale;
rad5=(rad+0.5)**2;
radp=rad+fExternRadius*fScale;
rad=rad+radprobe*fScale;
rad4=(rad+0.5)**2;
rad2=rad*rad;
radp2=radp*radp;

//set dielectric map
//check if sphere sits within limits of box
itest2=false;

//2011-05-13 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
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

if (itest2||itobig) { //slow method;
//2011-05-13 Seems redundant statement
rad2a = rad2 - 0.25;

for(iz=ismin.nZ;iz<=ismax.nZ;iz++){
for(iy=ismin.nY;iy<=ismax.nY;iy++){
for(ix=ismin.nX;ix<=ismax.nX;ix++){
//2011-05-13 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
ixyz=int_coord(ix,iy,iz);
dxyz=Int2Float(ixyz)-xn;
distsq=dxyz.dot.dxyz;
dxyz=dxyz+distsq;

if (dxyz.nX<rad2a) {
iEpsMap[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
}// if

if (dxyz.nY<rad2a) {
iEpsMap[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
}// if

if (dxyz.nZ<rad2a) {
iEpsMap[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
}// if

if(distsq<radp2) bDebMap[ix][iy][iz] =false;
}// do
}// do
}// do
else{ //faster method;
//IT HAS PROBLEMS!!!! Walter (be careful before using
//also in multidilectric case!!!&&!isitmd
rad2a=rad2-0.25;

//2011-05-13 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
ixn=nFloat2Int(xn);
fxn=Int2Float(ixn)-xn;
rad2av=rad2a-fxn;

for(ix=-lim;ix<=lim;ix++){
vtemp=Int2Float(ix)+fxn;
sq[ix]=vtemp*vtemp;
rad2aav[ix]=rad2a-vtemp;
}// do

//adjust inter-atom, different epsilon bgps+++04/2004 Walter
if (iNMedia>1&&bOnlyMol) {
for(i=1;i<=ibox;i++){
//2011-05-14 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
i123=ioff[i];
ixyz=ixn+i123;
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;
distsq=sq[i123.nX].nX+sq[i123.nY].nY+sq[i123.nZ].nZ;

if (distsq<rad2aav[i123.nX].nX) {
iac=mod(iEpsMap[ix][iy][iz].nX][epsdim)-1;

if (iac==-1||iac>iNatom) {
//occhio! non so cosa fa con i pori!!
iEpsMap[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
}
else{
//2011-05-14 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nX= ddxyz.nX+0.5;
dis2min1=ddxyz.dot.ddxyz-rad2;

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nX= ddxyz.nX+0.5;
dis2min2=ddxyz.dot.ddxyz-(sDelPhiPDB[iac].radius*fScale)**2;

if (dis2min2>dis2min1) iac=iv;
iEpsMap[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
}// if
}// if

if (distsq<rad2aav[i123.nY].nY) {
iac=mod(iEpsMap[ix][iy][iz].nY][epsdim)-1;
if (iac==-1||iac>iNatom) {
//occhio! non so cosa fa con i pori!!
iEpsMap[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
}
else{
//2011-05-14 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nY= ddxyz.nY+0.5;
dis2min1=ddxyz.dot.ddxyz-rad2;

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nY= ddxyz.nY+0.5;
dis2min2=ddxyz.dot.ddxyz-(sDelPhiPDB[iac].radius*fScale)**2;

if (dis2min2>dis2min1) iac=iv;
iEpsMap[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
}// if
}// if

if (distsq<rad2aav[i123.nZ].nZ) {
iac=mod(iEpsMap[ix][iy][iz].nZ][epsdim)-1;
if (iac==-1||iac>iNatom) {
//occhio! non so cosa fa con i pori!!
iEpsMap[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
}
else{
//2011-05-14 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
ddxyz=Int2Float(ixyz)-xn;
ddxyz.nZ= ddxyz.nZ+0.5;
dis2min1=ddxyz.dot.ddxyz-rad2;

ddxyz=Int2Float(ixyz)-xn2[iac];
ddxyz.nZ=ddxyz.nZ+0.5;
dis2min2=ddxyz.dot.ddxyz-(sDelPhiPDB[iac].radius*fScale)**2;

if (dis2min2>dis2min1) iac=iv;
iEpsMap[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
}// if
}// if

if(distsq<radp2) bDebMap[ix][iy][iz]=false;
}// do
}
else{
for(i=1;i<=ibox;i++){
//2011-05-14 Using operations on coord and int_coord
//type variables defined in module
//operators_on_coordinates
i123=ioff[i];
ixyz=ixn+i123;
ix=ixyz.nX;
iy=ixyz.nY;
iz=ixyz.nZ;
distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;

if (distsq<rad2aav[i123.nX].nX) {
iEpsMap[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
}// if

if (distsq<rad2aav[i123.nY].nY) {
iEpsMap[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
}// if

if (distsq<rad2aav[i123.nZ].nZ) {
iEpsMap[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
}// if

if (distsq<radp2) bDebMap[ix][iy][iz]=false;
}// do
}// if
}// if
}// do DoATOMS //end do of atoms;

//Array deallocation is made by ordinary F95 statement
if(allocated(ioff)) deallocate(ioff);

if (bDebug) {
open(52,file='iepsmapnewfirst',form='formatted');
write (52,*)iGrid,iGrid*iGrid;
for(ix=1;ix<=iGrid;ix++){
for(iy=1;iy<=iGrid;iy++){
iz=(iGrid+1)/2;
write (52][*)ix][iy][iz][iEpsMap[ix][iy][iz].nZ;
}// do
}// do
close (52);
}// if

if(bVerbose)cout <<"Ending creating Van der Waals Epsilon Map " << endl;

return;

}// void setout;

