      if (isolv) then
         !increased the starting dimension of crgatn and..+++++
         extracrg=0      
         if (ndistr.gt.0) extracrg=igrid**3
         !2011-05-30 Allocation of the arrays below is moved to the 
         !           body of crgarr subroutine, arrays are accessible 
         !           via pointers module. Sizes of arrays are 
         !           determined before allocation inside the crgarr 
         !           subroutine
         !i_crgatn=memalloc(i_crgatn,4,natom+extracrg)
         !i_chrgv2=memalloc(i_chrgv2,realsiz,4*(natom+extracrg))
         !i_nqgrdtonqass=memalloc(i_nqgrdtonqass,4,natom+extracrg)
         !i_atmeps=memalloc(i_atmeps,realsiz,natom+extracrg)
         !i_chgpos=memalloc(i_chgpos,realsiz,3*ncrgmx)
         !2011-05-30 Parameters transfered via qlog and pointers modules
         call crgarr
         if (logs.or.lognl) then
            ico=0
            do ic=1,nqass
               if ((chgpos(ic).vorlt.xl).or.(chgpos(ic).vorgt.xr)) then
                  if (crgatn(ic).lt.0) then
                     write (6,'(''!WARNING: distribution'',I4,&
                           &'' outside the box'')') (-crgatn(ic))
                  else
                     if (crgatn(ic).gt.natom) then
                        write(6,'(''WARNING:crg'',I4,'',object'',I4,&
                             &'' outside the box'',3f8.3)')&
                             &ic,(crgatn(ic)-natom),chgpos(ic)
                     else
                        write(6, '(''!!! WARNING : charge '',&
                             &a15,'' outside the box'')') &
                             &  delphipdb(crgatn(ic))%atinf
                     end if
                  end if
                  ico=1
               end if
            end do
            if (ico.gt.0.and.ibctyp.ne.3) then
               write(6,*)'CHARGES OUTSIDE THE BOX AND &
                         &NOT DOING FOCUSSING, THEREFORE STOP'
               stop
            end if
         end if
      end if

