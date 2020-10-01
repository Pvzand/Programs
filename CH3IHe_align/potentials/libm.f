c------------------------------ libm.f ---------------------------------
c
c      collection of subroutines extracted from Marius'lib to make a
c      portable ch3i-he_n package
c                                              m. lewerenz  11/feb/11
c
c-----------------------------------------------------------------------
c-------------------------- error handling -----------------------------

      subroutine errprt(iout,pgname,text,icode)
c
c      prints error messages from library subroutines
c
c      iout   : unit number for message output, 0-> no output; input
c      pgname : name of the subroutine calling errprt; input
c      text   : message text; input
c      icode  : severity code: 0 -> warning, < 0 -> fatal error with
c               abort, else -> error but execution continues
c      subroutines called : getstr                   m. lewerenz dec/93
c
      character pgname*(*),text*(*),header*20,tail*40
      save nerror,nwarn,icall
      common /errcnt/ maxerr,maxwrn
      data icall/0/
c
      if(icall.eq.0) then
        icall=1
        nerror=0
        nwarn=0
      end if
c
      if(icode.lt.0) then
        header='  *** fatal error,'
        tail=', execution aborted ***'
      else if(icode.eq.0) then
        header='  *** warning,'
        tail=' ***'
        nwarn=nwarn+1
      else
        header='  *** error,'
        tail=', return without action ***'
        nerror=nerror+1
      end if
c
c      write the message on unit iout
c
      if(iout.gt.0) then
        call getstr(pgname,lname,iout)
        call getstr(text,ltext,iout)
        call getstr(header,lhead,iout)
        call getstr(tail,ltail,iout)
        write(iout,'(/6a/)') header(1:lhead),' ',text(1:ltext),' in ',
     #                      pgname(1:lname),tail(1:ltail)
        call flush(iout)
      end if
c
      jcode=icode
      if(maxerr.gt.0.and.nerror.ge.maxerr) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of errors exceeded, program stopped *** '
        jcode=-1
      end if
      if(maxwrn.gt.0.and.nwarn.ge.maxwrn) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of warnings exceeded, program stopped ***'
        jcode=-1
      end if
      if(iout.gt.0) call flush(iout)
c
      if(jcode.lt.0) stop
      return
c
c      error report, returns current number of errors and warning 
c
      entry errnum(nerr,nwrn)
      nerr=nerror
      nwrn=nwarn
      return
      end


c----------------------------------------------------------------------


      subroutine getstr(string,ls,iout)
c
c      eliminates leading and trailing blanks from string and returns
c      length of remaining string
c
c      string : character string; input & output
c               on output leading and trailing blanks of the original
c               string are missing
c      ls     : length of non blank portion of string; output
c      iout   : fortran unit for messages, silent if iout.le.0; input
c      subroutines called: errprt                    m. lewerenz jul/00
c
      character string*(*)
c
      ls=len(string)
      if(ls.gt.0) then
        i=0
    5   i=i+1
        if(i.le.ls) then
          if(string(i:i).eq.' ') then
            goto 5
          else
            ls=ls+1-i
            if(i.gt.1) then
              do j=1,ls
                string(j:j)=string(j+i-1:j+i-1)
              end do
            end if
   10       continue
            if(string(ls:ls).eq.' ') then
              ls=ls-1
              goto 10
            end if
          end if
        else
          ls=0
          call errprt(iout,'getstr','empty string',0)
        end if
      end if
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine ivsum(iv,nv,isum)
c
c      sums all elements of integer vector iv in isum, vector
c      length nv.   m. lewerenz 14/jun/90
c
      dimension iv(nv)
      isum=0
      if(nv.le.0) return
      do 10 i=1,nv
      isum=isum+iv(i)
   10 continue
      return
      end
c
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c-----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine nextin(iunit,noblnk,items,iout)
c
c      positions unit iunit on the next input line
c      (lines starting with #,>,% are considered as comments)
c      and reports the number of input items on this line.
c
c      iunit  : fortran unit to be searched; input
c      noblnk : switch indicating treatment of blank lines:
c               noblnk=0 => blank lines are treated as input lines,
c               otherwise blank lines are also skipped; input
c      items  : number of blank or comma separated items found on the
c               next input line; output
c               items is -1 if end of file was found.
c               items is -2 if a reading error occurred
c      iout   : fortran unit for messages; silent for iout.le.0; input
c      subroutines called : errprt,strlen            m. lewerenz jul/00
c
      implicit real*8 (a-h,o-z)
      character line*256
c
      if(iunit.le.0) then
        call errprt(iout,'nextin','invalid fortran unit',1)
      else
   10   continue
        read(iunit,'(a)',end=998,err=999) line
        call strlen(line,ll,0)
        if(noblnk.ne.0.and.ll.le.0) goto 10
        if(line(1:1).eq.'#'.or.line(1:1).eq.'>'
     &                     .or.line(1:1).eq.'%') goto 10
c
        items=0
        icopy=0
        i=0
   15   i=i+1
        if(i.le.ll) then
          if(line(i:i).ne.' '.and.line(i:i).ne.',') then
            if(icopy.eq.0) then
              icopy=1
              items=items+1
            end if
          else
            icopy=0
          end if
          goto 15
        end if
        backspace(unit=iunit)
      end if
      return
c
  998 backspace(unit=iunit)
      call errprt(iout,'nextin','end of file reached',0)
      items=-1
      return
  999 call errprt(iout,'nextin','error reading file',1)
      items=-2
      return
      end
C
C----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
C
      SUBROUTINE STRLEN(STRING,LS,IOUT)
C
C      DETERMINES LENGTH OF STRING
C
      CHARACTER STRING*(*)
C
      LS=LEN(STRING)
   10 IF(LS.GT.0.AND.STRING(LS:LS).EQ.' ') THEN
        LS=LS-1
        GOTO 10
      END IF
      IF(LS.EQ.0) call errprt(iout,'strlen','empty string',0)
      RETURN
      END
c


