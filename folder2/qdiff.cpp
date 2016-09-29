#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;
if (isolv) {
//increased the starting dimension of crgatn and..+++++
extracrg=0;
if (ndistr>0) extracrg=iGrid**3;
//2011-05-30 Allocation of the arrays below is moved to the
// body of crgarr void, arrays are accessible
// via pointers module. Sizes of arrays are
// determined before allocation inside the crgarr
// void
//i_crgatn=memalloc(i_crgatn,4,iNatom+extracrg)
//i_chrgv2=memalloc(i_chrgv2,realsiz,4*(iNatom+extracrg))
//i_nqgrdtonqass=memalloc(i_nqgrdtonqass,4,iNatom+extracrg)
//i_atmeps=memalloc(i_atmeps,realsiz,iNatom+extracrg)
//i_chgpos=memalloc(i_chgpos,realsiz,3*ncrgmx)
//2011-05-30 Parameters transfered via qlog and pointers modules
call crgarr;
if (logs||lognl) {
ico=0;
for(ic=1;ic<=nqass;ic++){
if ((chgpos[ic].vorlt.xl)||(chgpos[ic].vorgt.xr)) {
if (crgatn[ic]<0) {
write (6]['(''//WARNING: distribution''][I4]['' outside the box'')') (-crgatn[ic]);
}
else{
if (crgatn[ic]>iNatom) {
write(6]['(''WARNING:crg''][I4][''][object''][I4][&'' outside the box''][3f8.3)')ic][(crgatn[ic]-iNatom)][chgpos[ic];
}
else{
write(6][ '(''//!! WARNING : charge ''][&a15]['' outside the box'')') sDelPhiPDB[crgatn[ic]]%atinf;
}// if
}// if
ico=1;
}// if
}// do
if (ico>0&&ibctyp!=3) {
cout <<"CHARGES OUTSIDE THE BOX AND NOT DOING FOCUSSING << THEREFORE STOP" << endl;
stop;
}// if
}// if
}// if

