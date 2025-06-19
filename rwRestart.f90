program rwRestart

  implicit none  

  real, parameter    :: var_fill = 999.

  type cmi2
     integer(4), pointer :: p(:,:)  => null() 
  end type cmi2
  type cmi3
     integer(4), pointer :: p(:,:,:)  => null() 
  end type cmi3
  type cmr1
     real(8), pointer :: p(:)  => null() 
  end type cmr1
  type cmr2
     real(8), pointer :: p(:,:)  => null() 
  end type cmr2
  type cmr3
     real(8), pointer :: p(:,:,:)  => null() 
  end type cmr3
  type cml1
     logical, pointer :: p(:)  => null() 
  end type cml1
  type cmi1
     integer(4), pointer :: p(:)  => null() 
  end type cmi1
  type cmr41
     real(4), pointer :: p(:)  => null() 
  end type cmr41
  type(cmi2),  pointer :: me1(:)
  type(cmi1),  pointer :: mcol(:),kh_fake(:)
  type(cmi3), pointer :: m3p(:)
  type(cmr3), pointer :: bndz(:)
  type(cmi2),  pointer :: ijk(:)
  type(cmr1), pointer :: ldp(:)!,hz(:)
  type(cmr2),  pointer :: pt3(:),pt2(:),pt3m(:),pt2m(:),pt3mm(:),pt2mm(:),bndice(:)!,iota(:)
  type(cmr41), pointer :: dp(:)
  type(cml1), pointer :: l21(:)

  character(11)  :: starta

  integer(4)     :: i,lun,ios ,narg
  integer(4)     :: ia
  integer(4)     :: t,n3d
  integer(4)     :: narea

  integer(4),     allocatable :: iwet2(:),iwet3(:)

  integer(4),     allocatable :: nlo(:),nla(:),km(:)
  real(8),        allocatable :: lo1(:),la1(:),dlo(:),dla(:)

  character(100), allocatable :: fnm_data(:),fnm_idx(:),fnm_depth(:)
  character(180)  :: filnam,restartnml

  integer, allocatable   :: ncid(:), time_id(:)
  integer, allocatable   :: tr3_id(:,:)
  integer, allocatable   :: tr2_id(:,:)
  logical :: ncfirst
  character(len=256), allocatable :: ncappend(:)
  character(len=256) :: ncfn
  integer :: nctimestep,cyear,iyea,cmonth,imon,cday,chour,iw2,iw3,j
  character(2) :: texhour
  character(2) :: chary
  real(8), allocatable :: lat(:),lon(:)
  real(8),allocatable :: nctim(:)
  
  integer(4)   :: iphig,iphim,ilamg,ilamm
  integer(4)   :: idphg,idphm,idlag,idlam
  integer(4)   :: ew, ns, layers, nzbnd, nubnd, nvbnd, nrivers,nza
  integer(4)   :: nudams, nvdams, nweirs
  real(8)      :: deltat,aush,rbot
  real(8)      :: rdphs,rdlas,rphis,rlams
  real(8)      :: latn,lonw,dlat,dlon

  real(4), allocatable :: dzi(:)

  namelist /dimensions/ ew, ns, layers, nzbnd, nubnd, nvbnd, nrivers,        &
                          nudams, nvdams, nweirs
  namelist /setup/      iphig,iphim,rphis,ilamg,ilamm,rlams,                 &
                          idphg,idphm,rdphs,idlag,idlam,rdlas,dzi,deltat,      &
                          aush,rbot

  character(len=256)      :: cfgfile_in, nidxfile_in, hzfile_in
  integer, dimension(9) :: nestinglevels_in, nestingto_in, timelevel_in
                          
  namelist /cfglist/ narea
  namelist /cfgfn/   cfgfile_in, nidxfile_in,  hzfile_in, timelevel_in, nestinglevels_in, nestingto_in 
  
  integer :: eu3,eu2
  integer,dimension(7) :: indx3
  integer,dimension(13) :: indx2,indx2_,indx2n,indx2h
  
  character(180)  :: restart_file,restart_out,mod_fname_1,mod_fname_2,nc_sname
  character(180),allocatable  :: mod_fname_3(:)
  logical :: create_nc,mod_restart,mod_arxive,iotamodel,read_restart
  character(6), dimension(7) :: names3d
  character(6), dimension(9) :: names2d
  integer, dimension(7) :: nc_3d,mod_3dparams
  integer, dimension(9) :: nc_2d,mod_2dparams
  integer, allocatable :: nc_domains(:),mod_domains(:)
  real :: mod_mindep,mod_change,sfill,tfill
  character(11)  :: time_default
  
  namelist /restart_options/ restart_file,restart_out,create_nc,names3d,nc_3d,names2d,nc_2d,nc_domains,&
             mod_fname_1,mod_fname_2,mod_fname_3,mod_restart,mod_arxive,mod_domains,mod_3dparams,mod_2dparams,&
             mod_mindep,mod_change,time_default,iotamodel,nc_sname,read_restart,sfill,tfill
             
  narg=command_argument_count()
  if(narg.gt.0)then
     call get_command_argument(1,restartnml)
  else
     print *,'No command line argument. Expecting initial rwRestart namelist file'
     stop 1       
  endif 
  lun = 1 
  open (unit=lun, file='cfg.nml')
  read (unit=lun, nml=cfglist)
  print *,'narea=',narea
  ! Allocate and get file names
  allocate(fnm_data(narea),fnm_idx(narea),fnm_depth(narea))
  allocate(nlo(narea),nla(narea),km(narea),iwet2(narea),iwet3(narea))
  allocate(ldp(narea))
  allocate(lo1(narea),la1(narea),dlo(narea),dla(narea))
  allocate(me1(narea),mcol(narea),kh_fake(narea),dp(narea)) !,hz(narea)
  do ia=1,narea
    cfgfile_in       = ' '
    nidxfile_in     = ' '
    hzfile_in     = ' '
    print *,' ia=',ia,' reading area ...'
    read (unit=lun, nml=cfgfn,IOSTAT=ios)
    print *,' ia=',ia,' sucess, ios=',ios
    if (ios/=0) call exit(9)
    fnm_data(ia)     = trim(cfgfile_in)
    fnm_idx(ia)     = trim(nidxfile_in)
    fnm_depth(ia)   = trim(hzfile_in)
    write(*,*) trim(fnm_data(ia)),' ',trim(fnm_idx(ia)),' ',trim(fnm_depth(ia))
    
    open(unit=2,file=trim(fnm_data(ia))//'.nml',status='OLD')
    read (2, nml=dimensions, iostat=ios)
    nla(ia)=ns
    nlo(ia)=ew
    km(ia)=layers
    if (ia==1) nza=nzbnd
    print *,ns,ew,layers
    allocate( dp(ia)%p(layers), ldp(ia)%p(layers),dzi(layers) )
    read (2, nml=setup, iostat=ios)
    close(2)
    
    latn=iphig+iphim/60.0+rphis/3600.
    lonw=ilamg+ilamm/60.0+rlams/3600.
    dlat=idphg+idphm/60.0+rdphs/3600.
    dlon=idlag+idlam/60.0+rdlas/3600.
    la1(ia)=latn
    lo1(ia)=lonw
    dla(ia)=dlat
    dlo(ia)=dlon
    dp(ia)%p(:)=dzi
    do i=1,layers
     if (i==1) then
      ldp(ia)%p(i)=0.5*dp(ia)%p(i)
     else
      ldp(ia)%p(i)=ldp(ia)%p(i-1)+0.5*(dp(ia)%p(i-1)+dp(ia)%p(i))
     end if
    end do
    !print *,ldp(ia)%p(1:6)
    deallocate(dzi)
    allocate(me1(ia)%p(0:nla(ia)+1,0:nlo(ia)+1))
    
    open(unit=2, file=trim(fnm_idx(ia)), form='unformatted', action='read',            &
         status='old', access='stream', iostat=ios, convert='big_endian')
    read(2,iostat=ios) me1(ia)%p(0:nla(ia)+1,0:nlo(ia)+1)
    !print *,'index is read',ios
    iw2 = 0
    do j = 1,ew
     do i = 1,ns
      if (me1(ia)%p(i,j)>0) iw2 = iw2+1
     enddo
    enddo
    allocate(mcol(ia)%p(0:iw2),kh_fake(ia)%p(0:iw2))
    read(2,iostat=ios) mcol(ia)%p(0:iw2)   
    !print *,'mcol read',ios 
    read(2,iostat=ios) kh_fake(ia)%p(0:iw2)
    !print *,'kh_fake read',ios
    iw3=0
    do i=1,iw2
     iw3=iw3+kh_fake(ia)%p(i) 
    end do 
    close(2,iostat=ios)
    print *,trim(fnm_idx(ia)),' file closed',ios,' iw2=',iw2,' iw3=',iw3
     
    iwet2(ia)=iw2
    iwet3(ia)=iw3
    
    !allocate(hz(ia)%p(0:iw3))  
    !open(2, file=trim(fnm_depth(ia)), form='unformatted', action='read',            &
    !     status='old', access='stream', iostat=ios, convert='big_endian')
    !read(2) hz(ia)%p(0:iw3)
    !close(2)
    
  enddo
  close(unit=lun)
  
  allocate(ncappend(narea),nc_domains(narea),mod_domains(narea),mod_fname_3(narea))
  
  open(unit=lun,file=restartnml)  !'rwRestart.nml'
  read (unit=lun, nml=restart_options)
  close(unit=lun)
  eu3=0
  do i=1,7
   if (nc_3d(i)>0) then
    print *,i,trim(names3d(i))
    eu3=eu3+1
    indx3(eu3)=i
   else if (mod_3dparams(i)>0) then
    print *,'3d parameter ',i,' from mod_3dparams with name ',names3d(i),' should be also included also in names3d. Exiting.'
    call exit(3)
   end if 
  end do
  eu2=0
  do i=1,9
   if (nc_2d(i)>0) then
    print *,i,trim(names2d(i))
    eu2=eu2+1
    indx2(eu2)=i
    indx2_(eu2)=i
    indx2n(eu2)=i
    indx2h(eu2)=1
    if (i>7) then
     indx2(eu2+1)=999
     indx2(eu2+2)=999
     indx2_(eu2)=8
     indx2_(eu2+1)=9
     indx2_(eu2+2)=10
     indx2n(eu2+1)=indx2n(eu2)
     indx2n(eu2+2)=indx2n(eu2)
     indx2h(eu2+1)=2
     indx2h(eu2+2)=3
     if (i>8) then
      indx2_(eu2)=11
      indx2_(eu2+1)=12
      indx2_(eu2+2)=13
     end if
     eu2=eu2+2
    end if
   else if (mod_2dparams(i)>0) then
    print *,'2d parameter ',i,' from mod_2dparams with name ',names2d(i),' should be also included also in names2d. Exiting.'
    call exit(2)
   end if 
  end do
  print *,eu3,eu2
  print *,indx3(1:eu3)
  print *,indx2(1:eu2)
  print *,indx2_(1:eu2)
  print *,indx2n(1:eu2)
  
  do ia=1,narea
     write(*,*)'nlo,nla,km=',nlo(ia),nla(ia),km(ia)
     write(chary,'(I2.2)') ia
     ncappend(ia)=trim(nc_sname)//chary//'_'    !trim(fnm_data(ia)(:))
  end do
  
  allocate(m3p(narea))
  allocate(ijk(narea))
  do ia=1,narea
     n3d=nla(ia)*nlo(ia)*km(ia)
     allocate(m3p(ia)%p(nla(ia),nlo(ia),km(ia)))
     allocate(ijk(ia)%p(3,1:n3d))
  end do
  deallocate(fnm_data,fnm_idx,fnm_depth)

  do ia=1,narea
     m3p(ia)%p(1:,1:,1:)=0
     ijk(ia)%p(1:,1:)=0
     call perm(nla(ia),nlo(ia),me1(ia)%p(1:nla(ia),1:nlo(ia)),kh_fake(ia)%p(1:iwet2(ia)), &
              m3p(ia)%p(1:,1:,1:),ijk(ia)%p(1:,1:))
  end do

  allocate(pt3(narea),pt2(narea),pt3m(narea),pt2m(narea),pt3mm(narea),pt2mm(narea),l21(narea))
  allocate(bndz(1),bndice(1)) !,iota(narea))
  do ia=1,narea
     allocate(pt3(ia)%p(1:7,0:iwet3(ia))) !3d
     allocate(pt2(ia)%p(1:13,0:iwet2(ia))) !2d
     allocate(pt3m(ia)%p(1:eu3,0:iwet3(ia))) !3d
     allocate(pt2m(ia)%p(1:eu2,0:iwet2(ia))) !2d
     allocate(pt3mm(ia)%p(1:eu3,0:iwet3(ia))) !3d
     allocate(pt2mm(ia)%p(1:eu2,0:iwet2(ia))) !2d
     allocate(l21(ia)%p(0:iwet2(ia)))
  end do
  allocate(bndz(1)%p(1:2,1:km(1),0:nza),bndice(1)%p(1:3,0:nza))

  do ia=1,narea
     pt3(ia)%p(1:,0:)=0.0_8
     pt2(ia)%p(1:,0:)=0.0_8
     pt3m(ia)%p(1:,0:)=0.0_8
     pt2m(ia)%p(1:,0:)=0.0_8
     pt3mm(ia)%p(1:,0:)=0.0_8
     pt2mm(ia)%p(1:,0:)=0.0_8
     l21(ia)%p(0:)=.false.
     pt3(ia)%p(3,0:)=sfill
     pt3(ia)%p(4,0:)=tfill
     do i=1,eu3
      if (indx3(i)==3) then
       pt3m(ia)%p(i,0:)=sfill
       pt3mm(ia)%p(i,0:)=sfill
      end if
      if (indx3(i)==4) then
       pt3m(ia)%p(i,0:)=tfill
       pt3mm(ia)%p(i,0:)=tfill
      end if
     end do
  end do
  starta=time_default
  bndz(1)%p(:,:,:)=0.0_8
  bndice(1)%p(:,:)=0.0_8

  !---Read Restart File---------------------------------------------------------
if (read_restart) then
  lun=20
  open(lun,file=trim(restart_file),form='unformatted',iostat=ios,                      &
       action='read',access='stream',status='old', CONVERT="BIG_ENDIAN")
  if (ios /= 0) then
     write(*,*) 'Could not open restart file'
     call exit(1)
  end if
  read(lun,iostat=ios) starta
  if (ios /= 0) then
     write(*,*) 'I/O problem - iostat is:',ios
     write(*,*) 'Could not read starta from aufpkt'
     call exit(2)
  endif
  write(*,*)' starta ',starta

  do ia=1,narea
        i=0
        do t=1,7
          read(lun,iostat=ios) pt3(ia)%p(t,0:iwet3(ia))
          if (nc_3d(t)>0) then
            i=i+1
            pt3m(ia)%p(i,0:iwet3(ia))=pt3(ia)%p(t,0:iwet3(ia))
            pt3mm(ia)%p(i,0:iwet3(ia))=pt3(ia)%p(t,0:iwet3(ia))
          end if
          if (ios.ne.0) then
              write(*,*)'error: could not read 3d field from restart file ',ia,t,names3d(t)
              call exit(3)
          end if
        end do
        i=0
        read(lun,iostat=ios) pt2(ia)%p(1,0:iwet2(ia))
        if (nc_2d(1)>0) then
         i=i+1
         pt2m(ia)%p(1,0:iwet2(ia))=pt2(ia)%p(1,0:iwet2(ia))
         pt2mm(ia)%p(1,0:iwet2(ia))=pt2(ia)%p(1,0:iwet2(ia))
        end if
        read(lun,iostat=ios) l21(ia)%p(0:iwet2(ia))
        do t=2,8
          if (t<8) then
           read(lun,iostat=ios) pt2(ia)%p(t,0:iwet2(ia))
           if (nc_2d(t)>0) then
            i=i+1
            pt2m(ia)%p(i,0:iwet2(ia))=pt2(ia)%p(t,0:iwet2(ia))
            pt2mm(ia)%p(i,0:iwet2(ia))=pt2(ia)%p(t,0:iwet2(ia))
           end if
          else if (t==8) then
           read(lun,iostat=ios) pt2(ia)%p(8:10,1:iwet2(ia))
           if (nc_2d(t)>0) then
            i=i+1
            pt2m(ia)%p(i:i+2,0:iwet2(ia))=pt2(ia)%p(t:t+2,0:iwet2(ia))
            pt2mm(ia)%p(i:i+2,0:iwet2(ia))=pt2(ia)%p(t:t+2,0:iwet2(ia))
            i=i+2
           end if
          end if      
          if (ios.ne.0) then
              write(*,*)'error: could not read 2d field from restart file ',ia,t,names2d(t)
              call exit(2)
          end if
        end do
        j=i

     !    copy from the read_restart subroutine of HBM 
     !    read(lun,iostat=ios)                                                 &
     !    u(ia)%p(0:iw3(ia)),v(ia)%p(0:iw3(ia)),s(ia)%p(0:iw3(ia)),            &
     !    t(ia)%p(0:iw3(ia)),tke(ia)%p(0:iw3(ia)),diss(ia)%p(0:iw3(ia)),       &
     !    avv(ia)%p(0:iw3(ia)),                                                &
     !    z(ia)%p(0:iw2(ia)),casus(ia)%p(0:iw2(ia)),                           &
     !    ice(ia)%p(1,0:iw2(ia)),ice(ia)%p(2,0:iw2(ia)),ice(ia)%p(3,0:iw2(ia)),&
     !    tsnei(ia)%p(0:iw2(ia)),ueis(ia)%p(0:iw2(ia)),veis(ia)%p(0:iw2(ia)),  &
     !    t_soil(ia)%p(1:nsoilpar,1:iw2(ia))

        write(*,*)'region:',ia
        do i=1,eu3
         write(*,*)'min/max('//trim(names3d(indx3(i)))//')=',minval(pt3(ia)%p(indx3(i),1:)),maxval(pt3(ia)%p(indx3(i),1:))
        end do
  enddo
  read(lun,iostat=ios) bndz(1)%p(1,1:km(1),0:),bndz(1)%p(2,1:km(1),0:)
  read(lun,iostat=ios) bndice(1)%p(1,0:), bndice(1)%p(2,0:),bndice(1)%p(3,0:)
  if (iotamodel) then
   do ia=1,narea
    read(lun,iostat=ios) pt2(ia)%p(11:13,0:iwet2(ia))
    i=j
    if (nc_2d(9)>0) then
            i=i+1
            pt2m(ia)%p(i:i+2,0:iwet2(ia))=pt2(ia)%p(11:13,0:iwet2(ia))
            pt2mm(ia)%p(i:i+2,0:iwet2(ia))=pt2(ia)%p(11:13,0:iwet2(ia)) 
    end if
   end do
  end if
  close(lun)
end if
  
  read(starta(2:5),'(I4)') cyear
  read(starta(6:7),'(I2)') cmonth
  read(starta(9:10),'(I2)') cday
  if (starta(11:11)=='0') texhour='00'
  if (starta(11:11)=='2') texhour='06'
  if (starta(11:11)=='4') texhour='12'
  if (starta(11:11)=='6') texhour='18'
  read(texhour,'(I2)') chour
  imon=mod((cmonth+9),12)
  iyea=cyear-imon/10
 
  if (mod_restart) then
    do ia=1,narea
      if (mod_domains(ia)>0) then
         allocate( lat(nla(ia)), lon(nlo(ia)))
         do i=1,nla(ia)
            lat(i)=la1(ia)-dble(i-1)*dla(ia)
         end do
         do i=1,nlo(ia)
            lon(i)=lo1(ia)+dble(i-1)*dlo(ia)
         end do
         filnam=trim(mod_fname_1)//starta(2:5)//trim(mod_fname_2)//starta(2:7)//trim(mod_fname_3(ia))
         if (.not. mod_arxive) filnam=trim(mod_fname_1)//trim(mod_fname_3(ia))
         print *,filnam         
         do i=1,eu3
          if (mod_3dparams(indx3(i))>0) then
           call read_cmems(filnam,nla(ia),nlo(ia),km(ia),lat,lon, &
                     ldp(ia)%p(1:),iwet3(ia),m3p(ia)%p(1:,1:,1:),names3d(indx3(i)),pt3m(ia)%p(i,1:))
           do j=1,iwet3(ia)
            if (pt3(ia)%p(indx3(i),j)<900 .and. pt3m(ia)%p(i,j)<900 .and. ldp(ia)%p(ijk(ia)%p(3,j))>mod_mindep) then
             pt3mm(ia)%p(i,j)=mod_change*pt3m(ia)%p(i,j)+(1-mod_change)*pt3(ia)%p(indx3(i),j)
            end if
           end do
          end if
         end do        
         
         do i=1,eu2
          if (mod_2dparams(indx2n(i))>0) then
           call read_2d(filnam,nla(ia),nlo(ia),indx2h(i),lat,lon, &
                     iwet2(ia),m3p(ia)%p(1:,1:,1),names2d(indx2n(i)),pt2m(ia)%p(i,1:))
           do j=1,iwet2(ia)
            if (pt2(ia)%p(indx2_(i),j)<900 .and. pt2m(ia)%p(i,j)<900) then
             pt2mm(ia)%p(i,j)=mod_change*pt2m(ia)%p(i,j)+(1-mod_change)*pt2(ia)%p(indx2_(i),j)
            end if
           end do
           if (indx2(i)==2 .or. indx2(i)==3 .or. indx2(i)==4) then
            do j=1,iwet2(ia)
             if (pt2mm(ia)%p(i,j)>0.001) l21(ia)%p(j)=.true.
            end do
           end if
          end if
         end do         
         
         deallocate(lat,lon)
      end if
    end do
     
    open(lun,file=trim(restart_out),form='unformatted',iostat=ios,                      &
       action='write',access='stream',status='replace', CONVERT="BIG_ENDIAN")
    if (ios /= 0) then
     write(*,*) 'Could not create restart file'
     call exit(1)
    end if
    write(lun,iostat=ios) starta
    if (ios /= 0) then
     write(*,*) 'I/O problem - iostat is:',ios
     write(*,*) 'Could not write starta to restart'
     call exit(2)
    endif
    do ia=1,narea
        i=0
        do t=1,7
          if (nc_3d(t)>0) then
            i=i+1
            write(lun,iostat=ios) pt3mm(ia)%p(i,0:iwet3(ia))
          else
            write(lun,iostat=ios) pt3(ia)%p(t,0:iwet3(ia))
          end if
          if (ios.ne.0) then
              write(*,*)'error: could not write 3d field in restart file ',ia,t,names3d(t)
              call exit(3)
          end if
        end do
        i=0
        if (nc_2d(1)>0) then
         i=i+1
         write(lun,iostat=ios) pt2mm(ia)%p(1,0:iwet2(ia))
        else
         write(lun,iostat=ios) pt2(ia)%p(1,0:iwet2(ia))
        end if
        write(lun,iostat=ios) l21(ia)%p(0:iwet2(ia))
        do t=2,8
          if (t<8) then       
           if (nc_2d(t)>0) then
            i=i+1
            write(lun,iostat=ios) pt2mm(ia)%p(i,0:iwet2(ia))
           else
            write(lun,iostat=ios) pt2(ia)%p(t,0:iwet2(ia))
           end if
          else if (t==8) then
           if (nc_2d(t)>0) then
            i=i+1
            write(lun,iostat=ios) pt2mm(ia)%p(i:i+2,1:iwet2(ia))
            i=i+2
           else
            write(lun,iostat=ios) pt2(ia)%p(t:t+2,1:iwet2(ia))
           end if          
          end if
          if (ios.ne.0) then
              write(*,*)'error: could not write 2d field to restart file ',ia,t,names2d(t)
              call exit(2)
          end if
        end do
        j=i

     !    copy from the read_restart subroutine of HBM 
     !    read(lun,iostat=ios)                                                 &
     !    u(ia)%p(0:iw3(ia)),v(ia)%p(0:iw3(ia)),s(ia)%p(0:iw3(ia)),            &
     !    t(ia)%p(0:iw3(ia)),tke(ia)%p(0:iw3(ia)),diss(ia)%p(0:iw3(ia)),       &
     !    avv(ia)%p(0:iw3(ia)),                                                &
     !    z(ia)%p(0:iw2(ia)),casus(ia)%p(0:iw2(ia)),                           &
     !    ice(ia)%p(1,0:iw2(ia)),ice(ia)%p(2,0:iw2(ia)),ice(ia)%p(3,0:iw2(ia)),&
     !    tsnei(ia)%p(0:iw2(ia)),ueis(ia)%p(0:iw2(ia)),veis(ia)%p(0:iw2(ia)),  &
     !    t_soil(ia)%p(1:nsoilpar,1:iw2(ia))

    enddo
    write(lun,iostat=ios) bndz(1)%p(1,1:km(1),0:),bndz(1)%p(2,1:km(1),0:)
    write(lun,iostat=ios) bndice(1)%p(1,0:), bndice(1)%p(2,0:),bndice(1)%p(3,0:)
    if (iotamodel) then
     do ia=1,narea
      if (nc_2d(9)>0) then
       i=j+1
       write(lun,iostat=ios) pt2mm(ia)%p(i:i+2,0:iwet2(ia))
      else
       write(lun,iostat=ios) pt2(ia)%p(11:13,0:iwet2(ia))    
      end if
     end do
    end if
    close(lun)     
  end if

if (create_nc) then

  allocate( ncid(narea), time_id(narea), tr3_id(narea,eu3), tr2_id(narea,eu2))
  allocate(nctim(1))
  ncfirst=.true.

  nctim(1)=cday + 365*iyea + iyea/4 - iyea/100 + iyea/400 + (imon*306 + 5)/10 - 693902
  nctim(1)=nctim(1)+chour/24.0_8
  print *,starta,cyear,cmonth,cday,chour,nctim(1)
  if (ncfirst) then
      do ia=1,narea
       if (nc_domains(ia)>0) then
          write(ncfn,'(a)') trim(ncappend(ia))//starta(2:7)//starta(9:10)//texhour//'.nc'
          call execute_command_line('rm -f '//trim(ncfn),exitstat=ios) 
          if (ios/=0) print *,'cannot remove ',ncfn,' ios=',ios
          write(*,*) 'Initializing file: '//trim(ncfn)
          allocate( lat(nla(ia)), lon(nlo(ia)))
          do i=1,nla(ia)
            lat(i)=la1(ia)-dble(i-1)*dla(ia)
          end do
          do i=1,nlo(ia)
            lon(i)=lo1(ia)+dble(i-1)*dlo(ia)
          end do
          call init_nc(ncfn,nla(ia),nlo(ia),km(ia),lat,lon, &
                     ldp(ia)%p(1:),ncid(ia),time_id(ia),starta(2:7)//starta(9:10)//texhour, &
                     tr3_id(ia,1:),tr2_id(ia,1:))
          deallocate(lat,lon)
          !print *,ncid(ia),time_id(ia),tr3_id(ia,1),tr3_id(ia,2),tr3_id(ia,3),tr2_id(ia,1),tr2_id(ia,2),tr2_id(ia,3)
       end if
      end do
      ncfirst=.false.
      nctimestep = 1
  end if
  write(*,*) 'Writing timestep no: ', nctimestep
  do ia=1,narea
      if (nc_domains(ia)>0) then
        write(*,*) '   area: ', ia
        call write_nc(ncid(ia),nla(ia),nlo(ia),km(ia),iwet3(ia),iwet2(ia),nctimestep, &
             time_id(ia),nctim,tr3_id(ia,1:),tr2_id(ia,1:),ijk(ia)%p(1:,1:),real(pt3(ia)%p(indx3(1:eu3),1:),4), &
               real(pt2(ia)%p(indx2_(1:eu2),1:),4))
      end if
  enddo
  nctimestep = nctimestep+1
     
  if (mod_restart) then
   nctim(1)=nctim(1)+1.0_8/24.0_8
   do ia=1,narea
      if (nc_domains(ia)>0) then
        call write_nc(ncid(ia),nla(ia),nlo(ia),km(ia),iwet3(ia),iwet2(ia),nctimestep, &
             time_id(ia),nctim,tr3_id(ia,1:),tr2_id(ia,1:),ijk(ia)%p(1:,1:),real(pt3m(ia)%p(1:,1:),4),real(pt2m(ia)%p(1:,1:),4)) 
      end if
   enddo  
   nctimestep = nctimestep+1
   nctim(1)=nctim(1)+1.0_8/24.0_8
   do ia=1,narea
      if (nc_domains(ia)>0) then
        call write_nc(ncid(ia),nla(ia),nlo(ia),km(ia),iwet3(ia),iwet2(ia),nctimestep, &
             time_id(ia),nctim,tr3_id(ia,1:),tr2_id(ia,1:),ijk(ia)%p(1:,1:),real(pt3mm(ia)%p(1:,1:),4),real(pt2mm(ia)%p(1:,1:),4)) 
      end if
   enddo  
  end if

  do ia=1,narea
    if (nc_domains(ia)>0) call close_nc(ncid(ia),starta(2:7)//starta(9:10))
  enddo
  deallocate(nctim)
  deallocate(ncid, time_id,tr3_id,tr2_id)
end if

  deallocate(nlo,nla,km,iwet2,iwet3)
  deallocate(lo1,la1,dlo,dla)
  deallocate(ncappend,nc_domains,mod_domains,mod_fname_3)
  do ia=1,narea
     deallocate(me1(ia)%p)
     deallocate(mcol(ia)%p)
     !deallocate(hz(ia)%p)
     deallocate(kh_fake(ia)%p)
     deallocate(m3p(ia)%p)
     deallocate(ijk(ia)%p)
     deallocate(dp(ia)%p)
     deallocate(ldp(ia)%p)
     deallocate(l21(ia)%p)
     deallocate(pt3(ia)%p)
     deallocate(pt3m(ia)%p)
     deallocate(pt3mm(ia)%p)
     deallocate(pt2(ia)%p)
     deallocate(pt2m(ia)%p)
     deallocate(pt2mm(ia)%p)
     !deallocate(iota(ia)%p)     
     !deallocate(pts(ia)%p)
  end do
  deallocate(bndz(1)%p,bndice(1)%p)
  deallocate(me1,mcol,m3p,ijk,pt3,pt2,pt3m,pt2m,pt3mm,pt2mm,l21,kh_fake,ldp,bndz,bndice) !,hz
  
contains
  
  subroutine perm(nla,nlo,me1,kh_fake,m3p,ijk)
  implicit none
   
    integer(4), intent(in)    :: nla,nlo
    integer(4), intent(in)    :: me1(1:,1:)
    integer(4), intent(in)    :: kh_fake(1:)
    integer(4), intent(inout) :: m3p(1:,1:,1:)
    integer(4), intent(inout) :: ijk(1:,1:)
    integer(4) :: i,j,k,mi,mis,kb

    mi=0
    do j=1,nlo
       do i=1,nla
          mis=me1(i,j)
          if (mis.gt.0) then
             mi=mi+1
             m3p(i,j,1)=mi
             ijk(1,mi)=i
             ijk(2,mi)=j
             ijk(3,mi)=1
          end if
       end do
    end do
    do j=1,nlo
     do i=1,nla
      mis=me1(i,j)
      if (mis>0) then
       kb=kh_fake(mis)
       do k=2,kb
                mi=mi+1
                m3p(i,j,k)=mi
                ijk(1,mi)=i
                ijk(2,mi)=j
                ijk(3,mi)=k
       end do
      end if
     end do
    end do 
  end subroutine perm

  subroutine init_nc(ncfn,nlat,nlon,ndep,lat,lon,dep, &
                     ncid,timeid,start_date, &
                     tr3,tr2)

    use netcdf, only : nf90_create, nf90_def_dim, nf90_def_var, nf90_put_att, &
                       nf90_enddef, nf90_put_var, nf90_clobber, NF90_UNLIMITED,&
                       NF90_REAL8, NF90_REAL, NF90_GLOBAL, nf90_netcdf4, &
                       nf90_def_var_deflate
    implicit none

    integer,            intent(in)  :: nlon, nlat, ndep
    character(len=256), intent(in)  :: ncfn
    character(len= 10), intent(in)  :: start_date
    real(8),            intent(in)  :: lat(nlat), lon(nlon), dep(ndep)
    integer,            intent(out) :: ncid
    integer,            intent(out) :: timeid
    integer,            intent(out) :: tr3(eu3),tr2(eu2)

    ! NetCDF attributes
    integer, parameter    :: defualt_cache_nelems = 1000
    integer, parameter    :: default_cache_size = 32000000
    real,    parameter    :: default_cache_preemption = 0.75

    ! Dimension attributes
    character (len = 256)  :: lat_name, lat_standard_name, lat_long_name
    character (len = 256)  :: lat_units, lat_axis
    character (len = 256)  :: lon_name, lon_standard_name, lon_long_name
    character (len = 256)  :: lon_units, lon_axis
    character (len = 256)  :: dep_name, dep3_name, dep_standard_name, dep_long_name
    character (len = 256)  :: dep_units, dep_axis, dep_unit_long, dep_positive
    character (len = 256)  :: time_name, time_standard_name, time_units
    character (len = 256)  :: time_calendar, time_long_name, time_axis
    character (len = 256)  :: time_CoordinateAxisType

    ! Variables for global attributes
    character (len = 256)  :: title, institution, references, contact, source
    character (len = 256)  :: Conventions, history, comment, grid_resolution

    ! Grid definition variables
    logical                :: nczip
    integer                :: tr
    integer, parameter     :: lun = 10
    integer                :: lon_dimid, lat_dimid, time_dimid
    integer                :: dep_dimid, dep3_dimid
    integer                :: lon_varid, lat_varid, dep_varid, dep3_varid
    integer                :: dimids4(4)
    character (len = 3)    :: ncformat
    real(8), dimension(3),parameter :: dep3=[1.,2.,3.]
    !real(8),allocatable    :: depe(:)

    namelist /nc_options/ncformat, nczip
    namelist /dim_atts/lat_name, lat_standard_name, &
                       lat_long_name, lat_units, lat_axis, &
                       lon_name, lon_standard_name, &
                       lon_long_name, lon_units, lon_axis, &
                       dep_name, dep_standard_name, &
                       dep_long_name, dep_units, dep_unit_long, dep_axis, &
                       dep_positive, &
                       time_name, time_standard_name, &
                       time_units, time_calendar, time_long_name, time_axis, &
                       time_CoordinateAxisType
    namelist /globatt/ title, institution, references, contact, source,        &
                       Conventions, history, comment, grid_resolution

    ! Namelist variables
    open(unit=lun,file='hbmnc.nml',status='old')

    ! Read nc format options
    read(unit=lun, nml=nc_options)

    ! Read dimension attributes
    read(unit=lun, nml=dim_atts)

    ! Read global attributes
    read(unit=lun, nml=globatt)
    close(lun)

    ! Create the ncdf file.
    if(ncformat=='nc') then
      call check( nf90_create(trim(ncfn), nf90_clobber, ncid) ) ! nc3 file
    elseif(ncformat=='nc4') then
      call check( nf90_create(trim(ncfn), nf90_netcdf4, ncid) ) ! nc4 file
    else
      write(*,*) 'ERROR!!! Invalid ncformat specified'
      return
    endif
  
   ! call check( nf90_create(trim(ncfn), nf90_netcdf4, ncid, &
   !               cache_size = default_cache_size,   & 
   !             cache_nelems = defualt_cache_nelems, &
   !         cache_preemption = default_cache_preemption))
    ! Define the dimensions.
    dep3_name='dep3'
    call check( nf90_def_dim(ncid, lat_name, nlat, lat_dimid) )
    call check( nf90_def_dim(ncid, lon_name, nlon, lon_dimid) )
    call check( nf90_def_dim(ncid, dep_name, ndep, dep_dimid) )
    call check( nf90_def_dim(ncid, dep3_name, 3, dep3_dimid) )
    call check( nf90_def_dim(ncid, time_name, NF90_UNLIMITED, time_dimid) )

    ! Define the coordinate variables.
    call check( nf90_def_var(ncid, lat_name, NF90_REAL8, lat_dimid, lat_varid) )
    call check( nf90_def_var(ncid, lon_name, NF90_REAL8, lon_dimid, lon_varid) )
    call check( nf90_def_var(ncid, dep_name, NF90_REAL8, dep_dimid, dep_varid) )
    call check( nf90_def_var(ncid, dep3_name, NF90_REAL8, dep3_dimid, dep3_varid) )
    call check( nf90_def_var(ncid, time_name,NF90_REAL8, time_dimid, timeid) )

    ! Compression settings
    if(nczip) then
      call check( nf90_def_var_deflate(ncid, lat_varid, shuffle = 1, &
                  deflate = 1, deflate_level = 3))
      call check( nf90_def_var_deflate(ncid, lon_varid, shuffle = 1, &
                  deflate = 1, deflate_level = 3))
      call check( nf90_def_var_deflate(ncid, dep_varid, shuffle = 1, &
                  deflate = 1, deflate_level = 3))
      call check( nf90_def_var_deflate(ncid, dep3_varid, shuffle = 1, &
                  deflate = 1, deflate_level = 3))
      call check( nf90_def_var_deflate(ncid, timeid, shuffle = 1, &
                  deflate = 1, deflate_level = 3))
    endif

    ! Assign units attributes to coordinate variables.
    call check( nf90_put_att(ncid, lat_varid, 'standard_name', &
                                           lat_standard_name))
    call check( nf90_put_att(ncid, lat_varid, 'long_name', lat_long_name) )
    call check( nf90_put_att(ncid, lat_varid, 'units', lat_units) )
    call check( nf90_put_att(ncid, lat_varid, 'axis', lat_axis) )
    call check( nf90_put_att(ncid, lon_varid, 'standard_name', &
                                           lon_standard_name))
    call check( nf90_put_att(ncid, lon_varid, 'long_name', lon_long_name) )
    call check( nf90_put_att(ncid, lon_varid, 'units', lon_units) )
    call check( nf90_put_att(ncid, lon_varid, 'axis', lon_axis) )
    call check( nf90_put_att(ncid, dep_varid, 'long_name', dep_long_name) )
    call check( nf90_put_att(ncid, dep_varid, 'units', dep_units) )
    call check( nf90_put_att(ncid, dep_varid, 'axis', dep_axis) )
    call check( nf90_put_att(ncid, dep3_varid, 'long_name', 'level') )
    call check( nf90_put_att(ncid, dep3_varid, 'units', '') )
    call check( nf90_put_att(ncid, dep3_varid, 'axis', 'z') )
    call check( nf90_put_att(ncid, timeid,'standard_name', time_standard_name))
    call check( nf90_put_att(ncid, timeid,'units', time_units) )
    call check( nf90_put_att(ncid, timeid,'calendar', time_calendar) )
    call check( nf90_put_att(ncid, timeid,'long_name', time_long_name) )
    call check( nf90_put_att(ncid, timeid,'axis', time_axis) )
    call check( nf90_put_att(ncid, timeid,'CoordinateAxisType',            &
                                          time_CoordinateAxisType) )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids4 = (/ lon_dimid, lat_dimid, dep_dimid, timeid /)
    ! Define the netCDF variables for the data.
!    call check( nf90_def_var(ncid, var_name, NF90_REAL, dimids, zid) )

    ! Assign units attributes to the netCDF variables.
    do tr=1,eu3
     call put_var_att(ncid,trim(names3d(indx3(tr))), 4, dimids4,  tr3(tr))
    end do

    do tr=1,eu2
     if (indx2(tr)<8) then
      call put_var_att(ncid,trim(names2d(indx2(tr))), 3, (/ lon_dimid, lat_dimid, timeid /),  tr2(tr))
     else if (indx2(tr)==8) then
      call put_var_att(ncid,trim(names2d(8)), 4, (/ lon_dimid, lat_dimid, dep3_dimid, timeid /),  tr2(tr))
     else if (indx2(tr)==9) then
      call put_var_att(ncid,trim(names2d(9)), 4, (/ lon_dimid, lat_dimid, dep3_dimid, timeid /),  tr2(tr))
     else
      tr2(tr)=-1
     end if
    end do

    ! Assign global attributes
    call check( nf90_put_att(ncid,NF90_GLOBAL,'title',title) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'institution',institution) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'references',references) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'contact',contact) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'source',source) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'Conventions',Conventions) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'history',history) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'comment',comment) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'grid_resolution', &
                                               grid_resolution) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'creation_date', 'Created '//    &
                                            currentdate()) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'start_date',start_date) )

    !print *,'3d ',trim(names3d(indx3(1))),' ',trim(names3d(indx3(2))),' ',trim(names3d(indx3(3))),' ',trim(names3d(indx3(eu3)))
    !print *,'2d ',trim(names2d(indx2(1))),' ',trim(names2d(indx2(2))),' ',trim(names2d(indx2(eu2)))
    ! End define mode.
    call check( nf90_enddef(ncid) )
    !print *,'aaa'

    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
    call check( nf90_put_var(ncid, lat_varid, lat) )
    !print *,'a1'
    call check( nf90_put_var(ncid, lon_varid, lon) )
    !print *,'a2'
    !print *,ncid,lat_dimid,lon_dimid,time_dimid,dep_dimid,lat_varid,lon_varid,timeid,dep_varid,ndep,&
    !         dep(1),dep(2),dep(ndep-1),dep(ndep)
    call check( nf90_put_var(ncid, dep_varid, dep) ) !dep_varid
    call check( nf90_put_var(ncid, dep3_varid, dep3) ) !dep_varid
    !print *,'bbb'
    
  end subroutine init_nc

! ---------------------------------------------------------------------------

  subroutine write_nc(ncid,nlat,nlon,ndep,iwet3,iwet2,timestep,time_varid,time, &
                      tr3,tr2,ij,valu3,valu2)

    use netcdf, only : nf90_put_var

    implicit none

    integer, intent(in)   :: nlon, nlat, ndep, iwet3, iwet2, timestep
    integer, intent(in)   :: ncid, time_varid
    real(8), intent(in)   :: time(1)
    integer, intent(in)   :: tr3(eu3),tr2(eu2)
    real, intent(in) :: valu3(eu3,1:iwet3)
    real, intent(in) :: valu2(eu2,1:iwet2)
    integer, intent(in)  :: ij(1:3,1:iwet3)
    real,dimension(:,:,:),allocatable :: work2
    real,dimension(:,:),allocatable :: work1
    integer :: i,tr

    call check( nf90_put_var(ncid,time_varid,time,start=(/ timestep /),count=(/ 1 /)) )

    ! Write 3D data.  
    allocate(work2(nlon,nlat,ndep))
    do tr=1,eu3
     work2(:,:,:)=var_fill
     !print *,tr,indx3(tr),valu3(1,indx3(tr)),valu3(iwet2,indx3(tr))!,valu3(iwet3,indx3(tr))
     do i=1,iwet3
      if (timestep==-1) then
       work2(ij(2,i),ij(1,i),ij(3,i))=valu3(indx3(tr),i)
      else
       work2(ij(2,i),ij(1,i),ij(3,i))=valu3(tr,i)
      end if
     end do
     !print *,'trying ',timestep,tr,tr3(tr),names3d(indx3(tr))
     call check( nf90_put_var(ncid,tr3(tr),work2(:,:,:),start=(/ 1, 1, 1, timestep /),count=(/ nlon, nlat, ndep, 1/)) )
     !print *,'success ',timestep,tr,tr3(tr),names3d(indx3(tr))
    end do
    deallocate(work2)

    ! Write 2D data.  
    allocate(work1(nlon,nlat))
    do tr=1,eu2
     if (indx2(tr)<8) then
     work1(:,:)=var_fill
     do i=1,iwet2
      if (timestep==-1) then
       work1(ij(2,i),ij(1,i))=valu2(indx2(tr),i)
      else
       work1(ij(2,i),ij(1,i))=valu2(tr,i)
      end if
     end do
     !print *,'trying ',timestep,tr,tr2(tr),names2d(indx2(tr))
     call check( nf90_put_var(ncid,tr2(tr),work1(:,:),start=(/ 1, 1, timestep /),count=(/ nlon, nlat, 1 /)) )
     !print *,'success ',timestep,tr,tr2(tr),names2d(indx2(tr))
     end if
    end do
    deallocate(work1)
    
    ! Write 2D data continued.  
    allocate(work2(nlon,nlat,3))
    do tr=1,eu2
     if (indx2(tr)==8 .or. indx2(tr)==9) then
      work2(:,:,:)=var_fill
      !print *,tr,indx3(tr),valu3(1,indx3(tr)),valu3(iwet2,indx3(tr))!,valu3(iwet3,indx3(tr))
      do i=1,iwet2
       if (timestep==-1) then
        work2(ij(2,i),ij(1,i),1)=valu2(indx2_(tr),i)
        work2(ij(2,i),ij(1,i),2)=valu2(indx2_(tr)+1,i)
        work2(ij(2,i),ij(1,i),3)=valu2(indx2_(tr)+2,i)
       else
        work2(ij(2,i),ij(1,i),1)=valu2(tr,i)
        work2(ij(2,i),ij(1,i),2)=valu2(tr+1,i)
        work2(ij(2,i),ij(1,i),3)=valu2(tr+2,i)
       end if
      end do
      !print *,'trying ',timestep,tr,tr3(tr),names3d(indx3(tr))
      call check( nf90_put_var(ncid,tr2(tr),work2(:,:,:),start=(/ 1, 1, 1, timestep /),count=(/ nlon, nlat, 3, 1/)) )
      !print *,'success ',timestep,tr,tr3(tr),names3d(indx3(tr))
     end if
    end do
    deallocate(work2)


   end subroutine write_nc
   
! ---------------------------------------------------------------------------

  subroutine put_var_att(ncid, varname, ndims, dimids, varid)

    use netcdf, only : nf90_def_var, nf90_put_att, NF90_REAL, &
                       nf90_def_var_deflate

    implicit none

    integer, intent(in)           :: ncid, ndims
    integer, intent(out)          :: varid
    integer                       :: dimids(ndims)
    character(len=*), intent(in)  :: varname

    ! Variable attributes
    logical              :: nczip
    character (len=3)    :: ncformat
    character (len=256)  :: var_units
    character (len=256)  :: units_long, coordinateaxes, standard_name, long_name
    real                 :: valid_min, valid_max

    namelist /nc_options/ncformat, nczip

    ! Namelist variables
    open(3,file='hbmnc.nml',status='old')

    ! Read nc format options
    read(3, nml=nc_options)
    close(3)

    ! Find variable attribute from namelist 
    standard_name=varname
    long_name=varname
    var_units=''
    units_long=''
    valid_min=-900
    valid_max=900
    var_units=''
    coordinateaxes = 'time_lon_lat' 

    ! Define the netCDF variable for the data.
    call check( nf90_def_var(ncid, varname, NF90_REAL, dimids, varid) )

    ! Compression settings
    if(nczip) then
      call check( nf90_def_var_deflate(ncid, varid, shuffle = 1, &
                  deflate = 1, deflate_level = 1))
    endif

    ! Assign units attributes to the netCDF variables.
    !call check( nf90_put_att(ncid, varid, '_CoordinateAxes', coordinateaxes) )
    call check( nf90_put_att(ncid, varid, 'standard_name', standard_name) )
    call check( nf90_put_att(ncid, varid, 'long_name', long_name) )
    call check( nf90_put_att(ncid, varid, 'units', var_units) )
    call check( nf90_put_att(ncid, varid, 'unit_long', units_long) )
    !call check( nf90_put_att(ncid, varid, 'valid_min', valid_min) )
    !call check( nf90_put_att(ncid, varid, 'valid_max', valid_max) )
    call check( nf90_put_att(ncid, varid, '_FillValue', var_fill) )
    call check( nf90_put_att(ncid, varid, 'missing_value', var_fill) )

 
  end subroutine put_var_att


  !-----------------------------------------------------------------------------
  
  subroutine read_cmems(filnam,nlat,nlon,ndep,lat,lon,ldp,iwet3,m3p,varname,valou)
  
  use netcdf, only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, &
                       nf90_inquire_dimension, nf90_nowrite

  implicit none
  
  integer, intent(in)   :: nlat,nlon,ndep,iwet3
  real(8), intent(in)  :: ldp(1:ndep)
  character(6), intent(in) :: varname
  character(180), intent(in) :: filnam
  real(8), intent(in) :: lat(nlat),lon(nlon)
  integer, intent(in)  :: m3p(nlat,nlon,ndep)
  real(8),intent(out)     :: valou(1:iwet3)
  integer :: ncid,ios,lat_varid,lon_varid,dep_varid,chl_varid,time_varid,nlats,nlons,ndeps,ntimes,i,j,k,xci,yci
  real(4), allocatable :: lats(:),lons(:),deps(:),chl_in(:,:,:)
  real(8)  :: sumsub,gribintbi,xxx,yyy
  integer, allocatable :: iklevel(:)
  real(8), allocatable :: works(:,:), workp(:,:),xcr(:),ycr(:),zzz(:)
  integer :: levas,levap,levs,levp
  logical :: istime
    
!    /net/isilon/ifs/arch/home/jis/CMEMS_MYP/2014/BALTICSEA_MULTIYEAR_PHY_003_011-TDS_201411.nc
  
  call check(nf90_open(filnam, nf90_nowrite, ncid))
  print *,'ncid=',ncid
  ios=nf90_inq_varid(ncid, 'lat', lat_varid)
  if (ios/=0) call check( nf90_inq_varid(ncid, 'latitude', lat_varid) )
  print *,'lat_varid',lat_varid
  ios=nf90_inq_varid(ncid, 'lon', lon_varid)
  if (ios/=0) call check( nf90_inq_varid(ncid, 'longitude', lon_varid) )
  print *,'lon_varid',lon_varid
  call check( nf90_inq_varid(ncid, 'depth', dep_varid) )
  print *,'dep_varid',dep_varid
  ios=nf90_inq_varid(ncid, 'time', time_varid)
  if (ios==0) then
   print *,'time_varid',time_varid
   istime=.true.
   call check( nf90_inquire_dimension(ncid, 1, len = ntimes) ) ! time_varid
   print *,ntimes
  else
   print *,'no time axis'
   istime=.false.
  end if
  print *,ia,varname
  call check( nf90_inq_varid(ncid, varname, chl_varid) )
  print *,'chl_varid',chl_varid
if (mod_arxive) then
  lat_varid=3
  lon_varid=4
  dep_varid=2
  !lat_varid=3
  !lon_varid=2
  !dep_varid=4
  call check( nf90_inquire_dimension(ncid, 3, len = nlats) ) !lat_varid
  print *,nlats
  call check( nf90_inquire_dimension(ncid, 4, len = nlons) ) !lon_varid
  print *,nlons
  call check( nf90_inquire_dimension(ncid, 2, len = ndeps) ) !dep_varid
  print *,ndeps
  allocate(lats(nlats),lons(nlons),deps(ndeps),chl_in(nlons,nlats,ndeps))
  print *,nlons,nlats,ndeps
  call check( nf90_get_var(ncid, 6, lats) )   !lat_varid
  print *,'ready1',lats(1),lats(2),lats(3),lats(nlats-2),lats(nlats-1),lats(nlats)
  call check( nf90_get_var(ncid, 3, lons) )   !lon_varid
  print *,'ready2',lons(1),lons(2),lons(3),lons(nlons-2),lons(nlons-1),lons(nlons)
  call check( nf90_get_var(ncid, 1, deps) )  !dep_varid
  print *,'ready3',deps(1),deps(2),deps(ndeps-1),deps(ndeps)
else
  call check( nf90_inquire_dimension(ncid, time_varid, len = ntimes) )
  call check( nf90_inquire_dimension(ncid, lat_varid, len = nlats) ) !lat_varid
  print *,nlats
  call check( nf90_inquire_dimension(ncid, lon_varid, len = nlons) ) !lon_varid
  print *,nlons
  call check( nf90_inquire_dimension(ncid, dep_varid, len = ndeps) ) !dep_varid
  print *,ndeps
  allocate(lats(nlats),lons(nlons),deps(ndeps),chl_in(nlons,nlats,ndeps))
  print *,nlons,nlats,ndeps
  call check( nf90_get_var(ncid, lat_varid, lats) )   !lat_varid
  print *,'ready1',lats(1),lats(2),lats(3),lats(nlats-2),lats(nlats-1),lats(nlats)
  call check( nf90_get_var(ncid, lon_varid, lons) )   !lon_varid
  print *,'ready2',lons(1),lons(2),lons(3),lons(nlons-2),lons(nlons-1),lons(nlons)
  call check( nf90_get_var(ncid, dep_varid, deps) )  !dep_varid
  print *,'ready3',deps(1),deps(2),deps(ndeps-1),deps(ndeps)
end if

  print *,ncid,lat_varid,lon_varid,dep_varid,chl_varid
  if (istime) then
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1, 1, ntimes /), &
                              count = (/ nlons, nlats, ndeps, ntimes /)) )
  else
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1, 1/), &
                              count = (/ nlons, nlats, ndeps/)) )  
  end if
  valou(:)=999.                            
  allocate(iklevel(ndep),zzz(ndep))
  do k=1,ndep
    i=1
    do while (ldp(k)>deps(i)+1e-3 .and. i<ndeps)
     i=i+1
    end do 
    iklevel(k)=i
    if (i==1) then
     zzz(k)=0
    else
     zzz(k)=(deps(i)-ldp(k))/(deps(i)-deps(i-1))
    end if
    if (zzz(k)<0) zzz(k)=0
    if (zzz(k)>1) zzz(k)=1
  end do
  allocate(xcr(nlat),ycr(nlon))
  do j = 1,nlon
    ycr(j)=(lon(j)-lons(1))*(nlons-1)/(lons(nlons)-lons(1))+1
  end do
  do i = 1,nlat
    xcr(i)=(lat(i)-lats(1))*(nlats-1)/(lats(nlats)-lats(1))+1
  end do
  allocate(works(nlons,nlats),workp(nlons,nlats))
  levs=-1
  levp=-1
  do k=1,ndep
    levas=iklevel(k)
    levap=levas
    if (iklevel(k)>1 .and. zzz(k)>1e-3) then
     levap=iklevel(k)-1
    end if
    if (levap/=levp) then
     if (levap==levs) then
      workp(:,:)=works(:,:)
     else
      do j=1,nlons
       do i=1,nlats
        workp(j,i)=chl_in(j,i,levap)
       end do
      end do
     end if
     levp=levap
    end if
    if (levas==levap .and. levs/=levp) then
     levs=levp
     works(:,:)=workp(:,:)
    else if (levas/=levs) then
     do j=1,nlons
      do i=1,nlats
       works(j,i)=chl_in(j,i,levas)
      end do
     end do
     levs=levas
    end if
    do j = 1,nlon
     yci=int(ycr(j))
     if (yci<1) yci=1
     if (yci>nlons-1) yci=nlons-1
     yyy=ycr(j)-yci
     if (yyy<0) yyy=0
     if (yyy>1) yyy=1
     do i = 1,nlat
      if (m3p(i,j,k)>0) then
      xci=int(xcr(i))
      if (xci<1) xci=1
      if (xci>nlats-1) xci=nlats-1
      xxx=xcr(i)-xci
      if (xxx<0) xxx=0
      if (xxx>1) xxx=1

      sumsub=0
      gribintbi=0
      if (abs(works(yci,xci))<100) then
       sumsub=sumsub+(1-yyy)*(1-xxx)*(1-zzz(k))
       gribintbi=gribintbi+(1-yyy)*(1-xxx)*(1-zzz(k))*works(yci,xci)
      end if
      if (abs(works(yci,xci+1))<100) then
       sumsub=sumsub+(1-yyy)*xxx*(1-zzz(k))
       gribintbi=gribintbi+(1-yyy)*xxx*(1-zzz(k))*works(yci,xci+1)
      end if
      if (abs(works(yci+1,xci))<100) then
       sumsub=sumsub+yyy*(1-xxx)*(1-zzz(k))
       gribintbi=gribintbi+yyy*(1-xxx)*(1-zzz(k))*works(yci+1,xci)
      end if
      if (abs(works(yci+1,xci+1))<100) then
       sumsub=sumsub+yyy*xxx*(1-zzz(k))
       gribintbi=gribintbi+yyy*xxx*(1-zzz(k))*works(yci+1,xci+1)
      end if
      if (abs(workp(yci,xci))<100) then
       sumsub=sumsub+(1-yyy)*(1-xxx)*zzz(k)
       gribintbi=gribintbi+(1-yyy)*(1-xxx)*zzz(k)*workp(yci,xci)
      end if
      if (abs(workp(yci,xci+1))<100) then
       sumsub=sumsub+(1-yyy)*xxx*zzz(k)
       gribintbi=gribintbi+(1-yyy)*xxx*zzz(k)*workp(yci,xci+1)
      end if
      if (abs(workp(yci+1,xci))<100) then
       sumsub=sumsub+yyy*(1-xxx)*zzz(k)
       gribintbi=gribintbi+yyy*(1-xxx)*zzz(k)*workp(yci+1,xci)
      end if
      if (abs(workp(yci+1,xci+1))<100) then
       sumsub=sumsub+yyy*xxx*zzz(k)
       gribintbi=gribintbi+yyy*xxx*zzz(k)*workp(yci+1,xci+1)
      end if
      if (sumsub<1e-6) then
       gribintbi=var_fill
      else
       gribintbi=gribintbi/sumsub
      end if
      if (abs(gribintbi)>100) then
       if (xci>1 .and. abs(workp(yci,xci-1))<100) then
        gribintbi=workp(yci,xci-1)
       else if (xci<nlats-1 .and. abs(workp(yci,xci+2))<100) then
        gribintbi=workp(yci,xci+2)
       else if (yci>1 .and. abs(workp(yci-1,xci))<100) then
        gribintbi=workp(yci-1,xci)
       else if (yci<nlons-1 .and. abs(workp(yci+2,xci))<100) then
        gribintbi=workp(yci+2,xci)
       else if (xci>1 .and. yci>1 .and. abs(workp(yci-1,xci-1))<100) then
        gribintbi=workp(yci-1,xci-1)
       else if (xci>1 .and. yci<nlons-1 .and. abs(workp(yci+2,xci-1))<100) then
        gribintbi=workp(yci+2,xci-1)
       else if (xci<nlats-1 .and. yci>1 .and. abs(workp(yci-1,xci+2))<100) then
        gribintbi=workp(yci-1,xci+2)
       else if (xci<nlats-1 .and. yci<nlons-1 .and. abs(workp(yci+2,xci+2))<100) then
        gribintbi=workp(yci+2,xci+2)
       else if (k>1 .and. m3p(i,j,k-1)>0 .and. abs(valou(m3p(i,j,k-1)))<100) then
        gribintbi=valou(m3p(i,j,k-1))
       else
        gribintbi=var_fill
       end if
      end if
!      print*,k,j,i,ycr(j),yci,yyy,xcr(i),xci,xxx,vali(i,j,k)
      valou(m3p(i,j,k))=gribintbi
     end if
     end do
    end do
  end do
  deallocate(chl_in)
  deallocate(works,workp,xcr,ycr,iklevel,zzz)
  call check( nf90_close(ncid) )

  end subroutine read_cmems

  !-----------------------------------------------------------------------------
  
  subroutine read_2d(filnam,nlat,nlon,ndep,lat,lon,iwet2,m3p,varname,valou)
  
  use netcdf, only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, &
                       nf90_inquire_dimension, nf90_nowrite

  implicit none
  
  integer, intent(in)   :: nlat,nlon,ndep,iwet2
  character(6), intent(in) :: varname
  character(180), intent(in) :: filnam
  real(8), intent(in) :: lat(nlat),lon(nlon)
  integer, intent(in)  :: m3p(nlat,nlon,1)
  real(8),intent(out)     :: valou(1:iwet2)
  integer :: ncid,ios,lat_varid,lon_varid,dep_varid,chl_varid,time_varid,nlats,nlons,ntimes,i,j,xci,yci
  real(4), allocatable :: lats(:),lons(:),chl_in(:,:,:)
  real(8)  :: sumsub,gribintbi,xxx,yyy
  real(8), allocatable :: works(:,:), xcr(:),ycr(:)
  logical :: istime,is3d
    
  call check(nf90_open(filnam, nf90_nowrite, ncid))
  print *,'ncid=',ncid
  ios=nf90_inq_varid(ncid, 'lat', lat_varid)
  if (ios/=0) call check( nf90_inq_varid(ncid, 'latitude', lat_varid) )
  print *,'lat_varid',lat_varid
  ios=nf90_inq_varid(ncid, 'lon', lon_varid)
  if (ios/=0) call check( nf90_inq_varid(ncid, 'longitude', lon_varid) )
  print *,'lon_varid',lon_varid
  ios=nf90_inq_varid(ncid, 'dep3', dep_varid)
  print *,'dep_varid',dep_varid,ios
  is3d=(ios==0)
  ios=nf90_inq_varid(ncid, 'time', time_varid)
  if (ios==0) then
   print *,'time_varid',time_varid
   istime=.true.
   call check( nf90_inquire_dimension(ncid, time_varid, len = ntimes) ) ! time_varid
   print *,ntimes
  else
   print *,'no time axis'
   istime=.false.
  end if
  print *,ia,varname
  call check( nf90_inq_varid(ncid, varname, chl_varid) )
  print *,'chl_varid',chl_varid
if (mod_arxive) then
  lat_varid=3
  lon_varid=4
  dep_varid=2
  call check( nf90_inquire_dimension(ncid, 3, len = nlats) ) !lat_varid
  print *,nlats
  call check( nf90_inquire_dimension(ncid, 4, len = nlons) ) !lon_varid
  print *,nlons
  allocate(lats(nlats),lons(nlons),chl_in(nlons,nlats,ndep))
  call check( nf90_get_var(ncid, 6, lats) )   !lat_varid
  print *,'ready1',lats(1),lats(2),lats(3),lats(nlats-2),lats(nlats-1),lats(nlats)
  call check( nf90_get_var(ncid, 3, lons) )   !lon_varid
  print *,'ready2',lons(1),lons(2),lons(3),lons(nlons-2),lons(nlons-1),lons(nlons)
else
  call check( nf90_inquire_dimension(ncid, time_varid, len = ntimes) )
  call check( nf90_inquire_dimension(ncid, lat_varid, len = nlats) ) !lat_varid
  print *,nlats
  call check( nf90_inquire_dimension(ncid, lon_varid, len = nlons) ) !lon_varid
  print *,nlons
  allocate(lats(nlats),lons(nlons),chl_in(nlons,nlats,ndep))
  call check( nf90_get_var(ncid, lat_varid, lats) )   !lat_varid
  print *,'ready1',lats(1),lats(2),lats(3),lats(nlats-2),lats(nlats-1),lats(nlats)
  call check( nf90_get_var(ncid, lon_varid, lons) )   !lon_varid
  print *,'ready2',lons(1),lons(2),lons(3),lons(nlons-2),lons(nlons-1),lons(nlons)
end if
  call check( nf90_inquire_dimension(ncid, dep_varid, len = i) ) !lon_varid
  print *,ncid,lat_varid,lon_varid,dep_varid,chl_varid,istime,is3d,ndep,ntimes,i
  if (istime .and. is3d) then
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1, 1, ntimes /), &
                              count = (/ nlons, nlats, ndep, ntimes /)) )
  else if (is3d) then
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1, 1/), &
                              count = (/ nlons, nlats, ndep/)) )  
  else if (istime) then
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1, ntimes/), &
                              count = (/ nlons, nlats, ntimes/)) )  
  else
   call check( nf90_get_var(ncid, chl_varid, chl_in, start = (/ 1, 1/), &
                              count = (/ nlons, nlats/)) ) 
  end if
  valou(:)=var_fill                        
  allocate(xcr(nlat),ycr(nlon))
  do j = 1,nlon
    ycr(j)=(lon(j)-lons(1))*(nlons-1)/(lons(nlons)-lons(1))+1
  end do
  do i = 1,nlat
    xcr(i)=(lat(i)-lats(1))*(nlats-1)/(lats(nlats)-lats(1))+1
  end do
  allocate(works(nlons,nlats))
     do j=1,nlons
      do i=1,nlats
       works(j,i)=chl_in(j,i,ndep)
      end do
     end do
    do j = 1,nlon
     yci=int(ycr(j))
     if (yci<1) yci=1
     if (yci>nlons-1) yci=nlons-1
     yyy=ycr(j)-yci
     if (yyy<0) yyy=0
     if (yyy>1) yyy=1
     do i = 1,nlat
      if (m3p(i,j,1)>0) then
      xci=int(xcr(i))
      if (xci<1) xci=1
      if (xci>nlats-1) xci=nlats-1
      xxx=xcr(i)-xci
      if (xxx<0) xxx=0
      if (xxx>1) xxx=1

      sumsub=0
      gribintbi=0
      if (abs(works(yci,xci))<100) then
       sumsub=sumsub+(1-yyy)*(1-xxx)
       gribintbi=gribintbi+(1-yyy)*(1-xxx)*works(yci,xci)
      end if
      if (abs(works(yci,xci+1))<100) then
       sumsub=sumsub+(1-yyy)*xxx
       gribintbi=gribintbi+(1-yyy)*xxx*works(yci,xci+1)
      end if
      if (abs(works(yci+1,xci))<100) then
       sumsub=sumsub+yyy*(1-xxx)
       gribintbi=gribintbi+yyy*(1-xxx)*works(yci+1,xci)
      end if
      if (abs(works(yci+1,xci+1))<100) then
       sumsub=sumsub+yyy*xxx
       gribintbi=gribintbi+yyy*xxx*works(yci+1,xci+1)
      end if
      if (sumsub<1e-6) then
       gribintbi=var_fill
      else
       gribintbi=gribintbi/sumsub
      end if
      if (abs(gribintbi)>100) then
       if (xci>1 .and. abs(works(yci,xci-1))<100) then
        gribintbi=works(yci,xci-1)
       else if (xci<nlats-1 .and. abs(works(yci,xci+2))<100) then
        gribintbi=works(yci,xci+2)
       else if (yci>1 .and. abs(works(yci-1,xci))<100) then
        gribintbi=works(yci-1,xci)
       else if (yci<nlons-1 .and. abs(works(yci+2,xci))<100) then
        gribintbi=works(yci+2,xci)
       else if (xci>1 .and. yci>1 .and. abs(works(yci-1,xci-1))<100) then
        gribintbi=works(yci-1,xci-1)
       else if (xci>1 .and. yci<nlons-1 .and. abs(works(yci+2,xci-1))<100) then
        gribintbi=works(yci+2,xci-1)
       else if (xci<nlats-1 .and. yci>1 .and. abs(works(yci-1,xci+2))<100) then
        gribintbi=works(yci-1,xci+2)
       else if (xci<nlats-1 .and. yci<nlons-1 .and. abs(works(yci+2,xci+2))<100) then
        gribintbi=works(yci+2,xci+2)
       else
        gribintbi=var_fill
       end if
      end if
      valou(m3p(i,j,1))=gribintbi
     end if
     end do
    end do
  deallocate(chl_in)
  deallocate(works,xcr,ycr)
  call check( nf90_close(ncid) )

  end subroutine read_2d    
  
  !----------------------------------------------------------------------------
  

  subroutine close_nc(ncid,end_date)

    use netcdf, only : nf90_close, nf90_redef, nf90_put_att, nf90_enddef, &
                       NF90_GLOBAL

    implicit none

    integer, intent(in)           :: ncid
    character(len=10), intent(in) :: end_date

    ! Add end date to global attributes
    call check( nf90_redef(ncid) )
    call check( nf90_put_att(ncid,NF90_GLOBAL,'end_date',end_date) )
    ! End define mode.
    call check( nf90_enddef(ncid) )

    ! Close the file. This causes netCDF to flush all
    ! buffers and make
    ! sure your data are really written to disk.
    call check( nf90_close(ncid) )

    write(*,*) 'NetCDF file closed !'

  end subroutine close_nc

  subroutine check(status)

    use netcdf, only : nf90_noerr, nf90_strerror

    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      write(*,*) trim(nf90_strerror(status))
      stop 'Stopped'
    end if
  end subroutine check

  character(29) function currentdate()

    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: val
    call date_and_time(date,time,zone,val)
    !call date_and_time(DATE=date,ZONE=zone)
    !call date_and_time(TIME=time)
    !call date_and_time(VALUES=val)

    write(currentdate,'(2(i2.2,"."),i4.4," ",2(i2.2,":"),i2.2," ",a5," UTC")') &
                                    val(3),val(2),val(1),val(5:6),val(7),zone

  end function currentdate


end program rwRestart
