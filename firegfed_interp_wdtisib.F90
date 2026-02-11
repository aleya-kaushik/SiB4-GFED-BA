!--------------------------------------------------------------
subroutine firegfed_interp(indx, lon, lat, &
           sibg)

!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv fire emissions
!     between their read times
!
!--------------------------------------------------------------

use kinds
use module_io, only: &
   fire_step, fire_seccur, fire_secnext
use module_oparams, only: &
    fire_leaff, fire_stwdf, &
    fire_cdbf, fire_metlf, fire_strlf
use module_param, only: &
    poolcon, firecon
use module_pparams, only: &
     mwc, mol_to_umol, &
     month_names, pdb, &
     drytoc
!use module_pftinfo, only: pft_num
use module_pftinfo
use module_poolinfo
use module_sib, only: gridcell_type
use module_sibconst, only: &
   npoolpft, npoollu, nsoil, ntpool, &
   fireb_print, fireb_stop, fireb_thresh
use module_time, only: &
   month, day, year, &
   dtisib, dtsib, sec_tot, wt_daily

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gridcell_type), intent(inout) :: sibg

!...parameters
integer(i4), parameter :: isave=5
real(r8) :: dnzero=1.E-10

!...interpolation variables
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: totemis, curemis, pcemis
real(r8) :: curemisc13, curemis_dm, curemisc13_dm

!...distribution variables
integer(i4) :: ntpft
integer(i4), dimension(:), allocatable :: tpref, tpnum
real(r4), dimension(:), allocatable :: tparea
!real(r8), dimension(:), allocatable :: tpagb, tpagbtemp
real(r8), dimension(:), allocatable :: tpbiomass
character(len=clen), dimension(:), allocatable :: tpname

real(r8) :: tfarea
integer(i4), dimension(:), allocatable :: tsortref
!real(r8), dimension(:), allocatable :: flossb
real(r8), dimension(:,:), allocatable :: flosspft, flosslu
real(r8), dimension(:,:), allocatable :: flosspft_dm, flosslu_dm
real(r8), dimension(:,:), allocatable :: fgainlu
!real(r8), dimension(:), allocatable :: flossbc13

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl
integer :: loc_c3g, loc_c4g
real(r8) :: myemis, myba, myfl, mymort, mytrans, tempemis, tempemisc13
real(r8) :: mort_fract, c_compl
real(r8) :: tmppooltot, tmppoolc13, tempemispoolc13, rcemispoolc13
real(r8) :: st_scalar, fl_scalar
integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
integer(i4) :: j,jref,k,kref,n,nref,m,mref
!real(r8) :: rcpooltest

!... scalars, more detail to be added later
st_scalar = 0.5 ! moisture stress scaler for mortality and combustion completeness
fl_scalar = 0 ! forest loss scaler for combustion completeness

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-6
frpc13 = pool_indx_froot_c13-6
crpc13 = pool_indx_croot_c13-6
wpc13 =  pool_indx_stwd_c13-6
ppc13 =  pool_indx_prod_c13-6

cdbp  = pool_indx_cdb-npoolpft/2
metlp = pool_indx_metl-npoolpft/2
strlp = pool_indx_strl-npoolpft/2
slitp = pool_indx_slit-npoolpft/2
slowp = pool_indx_slow-npoolpft/2
armp  = pool_indx_arm-npoolpft/2

cdbpc13  = pool_indx_cdb_c13 - npoolpft
metlpc13  = pool_indx_metl_c13 - npoolpft
strlpc13  = pool_indx_strl_c13 - npoolpft
slitpc13  = pool_indx_slit_c13 - npoolpft
slowpc13 = pool_indx_slow_c13 - npoolpft
armpc13  = pool_indx_arm_c13 - npoolpft


!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: fire_interp'
!print*,'poolpft_dloss(1,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(1,1),sibg%l(l)%poollt%poolpft_dloss(1,2),sibg%l(l)%poollt%poolpft_dloss(1,3)
!print*,'poolpft_dloss(2,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(2,1),sibg%l(l)%poollt%poolpft_dloss(2,2),sibg%l(l)%poollt%poolpft_dloss(2,3)
!print*,'poolpft(1) :', sibg%l(l)%poollt%poolpft(1)
!print*,'poolpft(2) :', sibg%l(l)%poollt%poolpft(2)
!print*,'poolpft(3) :', sibg%l(l)%poollt%poolpft(3)
!print*,' '
 


!-----------------------------------------------
! reset values
do l=1, sibg%g_nlu
   sibg%l(l)%poollt%resp_fire = dzero
   sibg%l(l)%poollt%resp_firec13 = dzero
   !sibg%l(l)%poollt%loss_fire_lay(1:5,:) = dzero
   !sibg%l(l)%pooldt%loss_fire_lay(1:6,:) = dzero
   sibg%l(l)%poollt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%poollt%loss_transf_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_transf_lay(:,:) = dzero
   sibg%l(l)%pooldt%gain_transf_lay(:,:) = dzero
enddo
   
! only continue if fire emissions are being used
IF (fire_step .le. izero) RETURN

! get scaling factors
facsibdrv = dble(fire_seccur-sec_tot) / dble(fire_step)

! only continue if fire emissions are valid at this time
IF ((facsibdrv .GT. 1) .OR. (facsibdrv .LT. 0)) THEN
   sibg%gprogt%ba_dbg = dzero
   sibg%gprogt%ba_en2 = dzero
   sibg%gprogt%ba_en3 = dzero
   sibg%gprogt%ba_dnf = dzero
   sibg%gprogt%ba_eb1 = dzero
   sibg%gprogt%ba_eb2 = dzero
   sibg%gprogt%ba_db1 = dzero
   sibg%gprogt%ba_db2 = dzero
   sibg%gprogt%ba_db3 = dzero
   sibg%gprogt%ba_shb = dzero
   sibg%gprogt%ba_sha = dzero
   sibg%gprogt%ba_c3a = dzero
   sibg%gprogt%ba_c3g = dzero
   sibg%gprogt%ba_c4g = dzero
   sibg%gprogt%ba_c3c = dzero
   sibg%gprogt%ba_c4c = dzero
   sibg%gprogt%ba_mze = dzero
   sibg%gprogt%ba_soy = dzero
   sibg%gprogt%ba_wwt = dzero
   RETURN
ENDIF

! interpolate BA
sibg%gprogt%ba_dbg = facsibdrv*sibg%gprogt%ba_dbg1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_dbg2
sibg%gprogt%ba_en2 = facsibdrv*sibg%gprogt%ba_en21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_en22
sibg%gprogt%ba_en3 = facsibdrv*sibg%gprogt%ba_en31 &
         + (1.-facsibdrv) * sibg%gprogt%ba_en32
sibg%gprogt%ba_dnf = facsibdrv*sibg%gprogt%ba_dnf1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_dnf2
sibg%gprogt%ba_eb1 = facsibdrv*sibg%gprogt%ba_eb11 &
         + (1.-facsibdrv) * sibg%gprogt%ba_eb12
sibg%gprogt%ba_eb2 = facsibdrv*sibg%gprogt%ba_eb21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_eb22
sibg%gprogt%ba_db1 = facsibdrv*sibg%gprogt%ba_db11 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db12
sibg%gprogt%ba_db2 = facsibdrv*sibg%gprogt%ba_db21 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db22
sibg%gprogt%ba_db3 = facsibdrv*sibg%gprogt%ba_db31 &
         + (1.-facsibdrv) * sibg%gprogt%ba_db32
sibg%gprogt%ba_shb = facsibdrv*sibg%gprogt%ba_shb1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_shb2
sibg%gprogt%ba_sha = facsibdrv*sibg%gprogt%ba_sha1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_sha2
sibg%gprogt%ba_c3a = facsibdrv*sibg%gprogt%ba_c3a1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3a2
sibg%gprogt%ba_c3g = facsibdrv*sibg%gprogt%ba_c3g1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3g2
sibg%gprogt%ba_c4g = facsibdrv*sibg%gprogt%ba_c4g1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c4g2
sibg%gprogt%ba_c3c = facsibdrv*sibg%gprogt%ba_c3c1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c3c2
sibg%gprogt%ba_c4c = facsibdrv*sibg%gprogt%ba_c4c1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_c4c2
sibg%gprogt%ba_mze = facsibdrv*sibg%gprogt%ba_mze1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_mze2
sibg%gprogt%ba_soy = facsibdrv*sibg%gprogt%ba_soy1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_soy2
sibg%gprogt%ba_wwt = facsibdrv*sibg%gprogt%ba_wwt1 &
         + (1.-facsibdrv) * sibg%gprogt%ba_wwt2

! new SiB-GFED fire calculation
! Emis = pool     * (MF*CC) * BA  * EF
! g    = kg DM m-2 *    %   * m2  * g(kgDM)^-1
! note: BA units are BA m2 PFT⁻¹ day⁻¹

! current code calculation:
! mol C m-2 = mol C m-2  *   fraction  * fraction
! BA units: BA m2 PFT⁻¹ day⁻¹: do we need a day->second conversion?
! tested this conversion in the firegfed_read_* files: results in values too low
 
   !totemis = sibg%gprogt%firec * dtsib
   ntpft = sibg%g_nlu
   ! number of land units/PFTs per cell -> =1 for Hyy when run as 1.0ENF
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft))
   !allocate(tpagb(ntpft),tpbiomass(ntpft))
   allocate(tpname(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft
   tpnum(:) = pft_num(tpref)
   do l=1,ntpft
     tpname(l) = pft_name(tpnum(l))
   enddo
   tparea(:) = sibg%l(1:ntpft)%larea

   !...remove carbon from all PFTs
   allocate(flosspft(ntpft,npoolpft))
   allocate(flosslu(ntpft,npoollu))
   allocate(fgainlu(ntpft,npoollu))
   allocate(flosspft_dm(ntpft,npoolpft))
   allocate(flosslu_dm(ntpft,npoollu))
   flosspft(:,:) = dzero
   flosslu(:,:) = dzero
   fgainlu(:,:) = dzero
   flosspft_dm(:,:) = dzero
   flosslu_dm(:,:) = dzero

   curemis = 0
   curemis_dm = 0
   curemisc13 = 0
   curemisc13_dm = 0

   do l=1,ntpft
      myl=l
      !... Iniital pool totals come from c13_iso_calc 
      !... for sum of (lp+wp+cbdp+metlp+strlp)
      tmppooltot = sibg%l(myl)%fract%poolemistotC
      tmppoolc13 = sibg%l(myl)%fract%poolemisc13

      !... read in BA based on pftname
      if (tpname(l) .eq. "dbg") myba=sibg%gprogt%ba_dbg
      if (tpname(l) .eq. "en2") myba=sibg%gprogt%ba_en2
      if (tpname(l) .eq. "en3") myba=sibg%gprogt%ba_en3
      if (tpname(l) .eq. "dnf") myba=sibg%gprogt%ba_dnf
      if (tpname(l) .eq. "eb1") myba=sibg%gprogt%ba_eb1
      if (tpname(l) .eq. "eb2") myba=sibg%gprogt%ba_eb2
      if (tpname(l) .eq. "db1") myba=sibg%gprogt%ba_db1
      if (tpname(l) .eq. "db2") myba=sibg%gprogt%ba_db2
      if (tpname(l) .eq. "db3") myba=sibg%gprogt%ba_db3
      if (tpname(l) .eq. "shb") myba=sibg%gprogt%ba_shb
      if (tpname(l) .eq. "sha") myba=sibg%gprogt%ba_sha
      if (tpname(l) .eq. "c3a") myba=sibg%gprogt%ba_c3a
      if (tpname(l) .eq. "c3g") myba=sibg%gprogt%ba_c3g
      if (tpname(l) .eq. "c4g") myba=sibg%gprogt%ba_c4g
      if (tpname(l) .eq. "c3c") myba=sibg%gprogt%ba_c3c
      if (tpname(l) .eq. "c4c") myba=sibg%gprogt%ba_c4c
      if (tpname(l) .eq. "mze") myba=sibg%gprogt%ba_mze
      if (tpname(l) .eq. "soy") myba=sibg%gprogt%ba_soy
      if (tpname(l) .eq. "wwt") myba=sibg%gprogt%ba_wwt

      !....remove C from leaf pool
      myfl = MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(lp) &
             - sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
             - poolcon(tpnum(myl))%poolpft_min(lp))
      mort_fract = firecon(myl)%mfl(lp) + &
                   st_scalar*(firecon(myl)%mfh(lp) - firecon(myl)%mfl(lp))
      c_compl = firecon(myl)%ccl(lp) + &
                   st_scalar*(firecon(myl)%cch(lp) - firecon(myl)%ccl(lp))
      mymort = myfl * myba * mort_fract
      myemis = mymort * c_compl
      mytrans = mymort * (1.0 - c_compl)
 
      sibg%l(myl)%poollt%loss_fire_lay(lp,1) = myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(lp,1) = mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
           + mymort
      flosspft(myl,lp) = myemis
      flosspft_dm(myl,lp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosspft(myl,lpc13) = sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis
      flosspft_dm(myl,lpc13) = sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis*drytoc
      sibg%l(myl)%poollt%loss_fire_lay(lpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(lpc13,1)*myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(lpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(lpc13,1)*mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(lpc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis*drytoc

      !....remove C from wood pool
      myfl = MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(wp) &
             - sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
             - poolcon(tpnum(myl))%poolpft_min(wp))
      mort_fract = firecon(myl)%mfl(wp) + &
                   st_scalar*(firecon(myl)%mfh(wp) - firecon(myl)%mfl(wp))
      c_compl = firecon(myl)%ccl(wp) + &
                   st_scalar*(firecon(myl)%cch(wp) - firecon(myl)%ccl(wp))
      mymort = myfl * myba * mort_fract
      myemis = mymort * c_compl
      mytrans = mymort * (1.0 - c_compl)

      sibg%l(myl)%poollt%loss_fire_lay(wp,1) = myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(wp,1) = mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
           + mymort
      flosspft(myl,wp) = myemis
      flosspft_dm(myl,wp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosspft(myl,wpc13) = sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis
      flosspft_dm(myl,wpc13) = sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis*drytoc
      sibg%l(myl)%poollt%loss_fire_lay(wpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(wpc13,1)*myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(wpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(wpc13,1)*mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(wpc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis*drytoc

      !....remove C emis from metabolic litter, gain from fire xfr
      myfl = MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(metlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(metlp,1))
      !mort_fract = firecon(myl)%mfl(metlp) + &
      !             st_scalar*(firecon(myl)%mfh(metlp) - firecon(myl)%mfl(metlp))
      c_compl = firecon(myl)%ccl(metlp) + &
                   st_scalar*(firecon(myl)%cch(metlp) - firecon(myl)%ccl(metlp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      sibg%l(myl)%pooldt%loss_fire_lay(metlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(metlp,1) = mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlp,1) &
           + mymort
      flosslu(myl,metlp) = myemis
      flosslu_dm(myl,metlp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosslu(myl,metlpc13) = sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis
      flosslu_dm(myl,metlpc13) = &
                    sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis*drytoc
      sibg%l(myl)%pooldt%loss_fire_lay(metlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(metlpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(metlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(metlpc13,1)*mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(metlpc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis*drytoc

      !....remove C from structural litter, gain from fire xfr
      myfl = MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(strlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(strlp,1))
      !mort_fract = firecon(myl)%mfl(strlp) + &
      !             st_scalar*(firecon(myl)%mfh(strlp) - firecon(myl)%mfl(strlp))
      c_compl = firecon(myl)%ccl(strlp) + &
                   st_scalar*(firecon(myl)%cch(strlp) - firecon(myl)%ccl(strlp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      sibg%l(myl)%pooldt%loss_fire_lay(strlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(strlp,1) = mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlp,1) &
           + mymort
      flosslu(myl,strlp) = myemis
      flosslu_dm(myl,strlp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosslu(myl,strlpc13) = sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis
      flosslu_dm(myl,strlpc13) = &
                    sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis*drytoc
      sibg%l(myl)%pooldt%loss_fire_lay(strlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(strlpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(strlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(strlpc13,1)*mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(strlpc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis*drytoc

      !....remove C from coarse dead biomass, gain from fire xfr
      myfl = MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(cdbp) &
           - sibg%l(myl)%pooldt%poollu_dloss(cdbp,1))
      !mort_fract = firecon(myl)%mfl(cdbp) + &
      !             st_scalar*(firecon(myl)%mfh(cdbp) - firecon(myl)%mfl(cdbp))
      c_compl = firecon(myl)%ccl(cdbp) + &
                   st_scalar*(firecon(myl)%cch(cdbp) - firecon(myl)%ccl(cdbp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      sibg%l(myl)%pooldt%loss_fire_lay(cdbp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(cdbp,1) = mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) &
           + mymort
      flosslu(myl,cdbp) = myemis
      flosslu_dm(myl,cdbp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosslu(myl,cdbpc13) = sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis
      flosslu_dm(myl,cdbpc13) = &
                    sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis*drytoc
      sibg%l(myl)%pooldt%loss_fire_lay(cdbpc13,1) = &
             sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%loss_transf_lay(cdbpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc13,1)*mytrans * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis*drytoc

      !...remove C from product
      myfl = MAX(dzero, &
              sibg%l(myl)%poollt%poolpft(pp) &
              - sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
              - poolcon(tpnum(myl))%poolpft_min(pp))
      mort_fract = firecon(myl)%mfl(pp) + &
                   st_scalar*(firecon(myl)%mfh(pp) - firecon(myl)%mfl(pp))
      c_compl = firecon(myl)%ccl(pp) + &
                   st_scalar*(firecon(myl)%cch(pp) - firecon(myl)%ccl(pp))
      mymort = myfl * myba * mort_fract
      myemis = mymort * c_compl
      mytrans = mymort * (1.0 - c_compl)

      sibg%l(myl)%poollt%loss_fire_lay(pp,1) = myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(pp,1) = mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(pp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
           + mymort
      flosspft(myl,pp) = myemis
      flosspft_dm(myl,pp) = myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc

      flosspft(myl,ppc13) = sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis
      flosspft_dm(myl,ppc13) = &
                    sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis*drytoc
      sibg%l(myl)%poollt%loss_fire_lay(ppc13,1) = &
            sibg%l(myl)%poollt%rcpoolpft_lay(ppc13,1)*myemis * dtisib
      sibg%l(myl)%poollt%loss_transf_lay(ppc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(ppc13,1)*mytrans * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) = &
          sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) &
          + sibg%l(myl)%poollt%rcpoolpft_lay(ppc13,1)*mymort
      curemisc13 = curemisc13 + sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis*drytoc
      
      !...remove C from roots
      myfl = MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(frp) &
           - sum(sibg%l(myl)%poollt%poolpft_dloss(frp,:)) &
           - poolcon(tpnum(myl))%poolpft_min(frp))
      mort_fract = firecon(myl)%mfl(frp) + &
                   st_scalar*(firecon(myl)%mfh(frp) - firecon(myl)%mfl(frp))
      c_compl = firecon(myl)%ccl(frp) + &
                   st_scalar*(firecon(myl)%cch(frp) - firecon(myl)%ccl(frp))
      mymort = myfl * myba * mort_fract
      myemis = mymort * c_compl
      mytrans = mymort * (1.0 - c_compl)

      DO s=1,nsoil 
         sibg%l(myl)%poollt%loss_fire_lay(frp,s) = &
              myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(frp,s)
         sibg%l(myl)%poollt%loss_transf_lay(frp,s) = &
              mytrans * dtisib * sibg%l(myl)%poollt%poolpft_flay(frp,s)
         sibg%l(myl)%poollt%poolpft_dloss(frp,s) = &
              sibg%l(myl)%poollt%poolpft_dloss(frp,s) &
              + mymort * sibg%l(myl)%poollt%poolpft_flay(frp,s)

         sibg%l(myl)%poollt%loss_fire_lay(frpc13,s) = &
              sibg%l(myl)%poollt%rcpoolpft_lay(frpc13,s)*myemis &
              * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
         sibg%l(myl)%poollt%loss_transf_lay(frpc13,s) = &
              sibg%l(myl)%poollt%rcpoolpft_lay(frpc13,s)*mytrans &
              * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
         sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) = &
              sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) &
              + sibg%l(myl)%poollt%rcpoolpft_lay(frpc13,s)*mymort &
              * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
      ENDDO
      flosspft(myl,frp) = myemis
      flosspft_dm(myl,frp) = myemis*drytoc
      flosspft(myl,frpc13) = sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis
      flosspft_dm(myl,frpc13) = &
                    sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc
      curemisc13 = curemisc13 + sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis*drytoc
           
      myfl = MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(crp) &
           - sum(sibg%l(myl)%poollt%poolpft_dloss(crp,:)) & 
           - poolcon(tpnum(myl))%poolpft_min(crp))
      mort_fract = firecon(myl)%mfl(crp) + &
                   st_scalar*(firecon(myl)%mfh(crp) - firecon(myl)%mfl(crp))
      c_compl = firecon(myl)%ccl(crp) + &
                   st_scalar*(firecon(myl)%cch(crp) - firecon(myl)%ccl(crp))
      mymort = myfl * myba * mort_fract
      myemis = mymort * c_compl
      mytrans = mymort * (1.0 - c_compl)

      DO s=1,nsoil 
         sibg%l(myl)%poollt%loss_fire_lay(crp,s) = &
              myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(crp,s)
         sibg%l(myl)%poollt%loss_transf_lay(crp,s) = &
              mytrans * dtisib * sibg%l(myl)%poollt%poolpft_flay(crp,s)
         sibg%l(myl)%poollt%poolpft_dloss(crp,s) = &
              sibg%l(myl)%poollt%poolpft_dloss(crp,s) &
              + mymort * sibg%l(myl)%poollt%poolpft_flay(crp,s)

         sibg%l(myl)%poollt%loss_fire_lay(crpc13,s) = &
              sibg%l(myl)%poollt%rcpoolpft_lay(crpc13,s)*myemis &
              * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
         sibg%l(myl)%poollt%loss_transf_lay(crpc13,s) = &
              sibg%l(myl)%poollt%rcpoolpft_lay(crpc13,s)*mytrans &
              * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
         sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) = &
              sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) &
              + sibg%l(myl)%poollt%rcpoolpft_lay(crpc13,s)*mymort &
              * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
      ENDDO
      flosspft(myl,crp) = myemis
      flosspft_dm(myl,crp) = myemis*drytoc
      flosspft(myl,crpc13) = sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis
      flosspft_dm(myl,crpc13) = &
                    sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc
      curemisc13 = curemisc13 + sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis*drytoc

      !...remove C from soil litter
      myfl = MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(slitp) &
           - sum(sibg%l(myl)%pooldt%poollu_dloss(slitp,:)))
      !mort_fract = firecon(myl)%mfl(slitp) + &
      !             st_scalar*(firecon(myl)%mfh(slitp) - firecon(myl)%mfl(slitp))
      c_compl = firecon(myl)%ccl(slitp) + &
                   st_scalar*(firecon(myl)%cch(slitp) - firecon(myl)%ccl(slitp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      DO s=1,nsoil 
        sibg%l(myl)%pooldt%loss_fire_lay(slitp,s) = &
             myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitp,s)
        sibg%l(myl)%pooldt%loss_transf_lay(slitp,s) = &
             mytrans * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitp,s)
        sibg%l(myl)%pooldt%poollu_dloss(slitp,s) = &
             sibg%l(myl)%pooldt%poollu_dloss(slitp,s) &
             + mymort * sibg%l(myl)%pooldt%poollu_flay(slitp,s)

        sibg%l(myl)%pooldt%loss_fire_lay(slitpc13,s) = &
             sibg%l(myl)%pooldt%rcpoollu_lay(slitpc13,s)*myemis &
             * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
        sibg%l(myl)%pooldt%loss_transf_lay(slitpc13,s) = &
              sibg%l(myl)%pooldt%rcpoollu_lay(slitpc13,s)*mytrans &
              * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
        sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) = &
             sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) &
             + sibg%l(myl)%pooldt%rcpoollu_lay(slitpc13,s)*mymort &
             * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
      ENDDO
      flosslu(myl,slitp) = myemis
      flosslu_dm(myl,slitp) = myemis*drytoc
      flosslu(myl,slitpc13) = sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis
      flosslu_dm(myl,slitpc13) = &
                  sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis*drytoc

      !...remove C from soil slow
      myfl = MAX(dzero, &
          sibg%l(myl)%pooldt%poollu(slowp) &
          - sum(sibg%l(myl)%pooldt%poollu_dloss(slowp,:)))
      !mort_fract = firecon(myl)%mfl(slowp) + &
      !             st_scalar*(firecon(myl)%mfh(slowp) - firecon(myl)%mfl(slowp))
      c_compl = firecon(myl)%ccl(slowp) + &
                   st_scalar*(firecon(myl)%cch(slowp) - firecon(myl)%ccl(slowp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      DO s=1,nsoil 
         sibg%l(myl)%pooldt%loss_fire_lay(slowp,s) = &
              myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowp,s)
         sibg%l(myl)%pooldt%loss_transf_lay(slowp,s) = &
              mytrans * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowp,s)
         sibg%l(myl)%pooldt%poollu_dloss(slowp,s) = &
              sibg%l(myl)%pooldt%poollu_dloss(slowp,s) &
              + mymort * sibg%l(myl)%pooldt%poollu_flay(slowp,s)

         sibg%l(myl)%pooldt%loss_fire_lay(slowpc13,s) = &
              sibg%l(myl)%pooldt%rcpoollu_lay(slowpc13,s)*myemis &
              * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
         sibg%l(myl)%pooldt%loss_transf_lay(slowpc13,s) = &
              sibg%l(myl)%pooldt%rcpoollu_lay(slowpc13,s)*mytrans &
              * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
         sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) = &
              sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) &
              + sibg%l(myl)%pooldt%rcpoollu_lay(slowpc13,s)*mymort &
              * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
      ENDDO
      flosslu(myl,slowp) = myemis
      flosslu_dm(myl,slowp) = myemis*drytoc
      flosslu(myl,slowpc13) = sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis
      flosslu_dm(myl,slowpc13) = &
                  sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis*drytoc

      !...remove C from soil passive(armored)
      myfl = MIN(curemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(armp) &
           - sum(sibg%l(myl)%pooldt%poollu_dloss(armp,:))))
      !mort_fract = firecon(myl)%mfl(armp) + &
      !             st_scalar*(firecon(myl)%mfh(armp) - firecon(myl)%mfl(armp))
      c_compl = firecon(myl)%ccl(armp) + &
                   st_scalar*(firecon(myl)%cch(armp) - firecon(myl)%ccl(armp))
      mymort = 0.0
      myemis = myfl * myba * c_compl
      mytrans = 0.0

      DO s=1,nsoil 
          sibg%l(myl)%pooldt%loss_fire_lay(armp,s) = &
               myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(armp,s)
         sibg%l(myl)%pooldt%loss_transf_lay(armp,s) = &
              mytrans * dtisib * sibg%l(myl)%pooldt%poollu_flay(armp,s)
          sibg%l(myl)%pooldt%poollu_dloss(armp,s) = &
               sibg%l(myl)%pooldt%poollu_dloss(armp,s) &
               + mymort * sibg%l(myl)%pooldt%poollu_flay(armp,s)

          sibg%l(myl)%pooldt%loss_fire_lay(armpc13,s) = &
               sibg%l(myl)%pooldt%rcpoollu_lay(armpc13,s)*myemis &
               * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
         sibg%l(myl)%pooldt%loss_transf_lay(armpc13,s) = &
              sibg%l(myl)%pooldt%rcpoollu_lay(armpc13,s)*mytrans &
              * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
          sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) = &
               sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) &
               + sibg%l(myl)%pooldt%rcpoollu_lay(armpc13,s)*mymort &
               * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
      ENDDO
      flosslu(myl,armp) = myemis
      flosslu_dm(myl,armp) = myemis*drytoc
      flosslu(myl,armpc13) = sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis
      flosslu_dm(myl,armpc13) = &
                  sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis*drytoc
      curemis = curemis + myemis
      curemis_dm = curemis_dm + myemis*drytoc
      curemisc13 = curemisc13 + sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis
      curemisc13_dm = curemisc13_dm + &
                      sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis*drytoc

      !-----Carbon Transfers from live pools-----
      !...n is the sending/from pool
      !...m is the receiving/to pool
      do n=1,npoolpft/2 !npoolpft was 5, now 10 with C-13 pools, (1,5) live total C
         !do m=npoolpft+1,ntpool !was (6,11) for dead pools 
         do m=npoolpft/2+1,ntpool/2 !now (6,11) 
            mref=m-npoolpft/2 !(m-5), i.e. 1-6 for npoollu dead pools
              if (poolcon(myl)%pool_trans_frac(n,m) > dzero) then
                !.....transfer from single-layer canopy/soil to surface
                if ((pool_indx_lay(n) .eq. 1) .and. &
                    (pool_indx_lay(m) .eq. 1)) then
                       sibg%l(myl)%pooldt%gain_transf_lay(mref,1) = &
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,1) + &
                           sibg%l(myl)%poollt%loss_transf_lay(n,1) &
                           * poolcon(myl)%pool_trans_frac(n,m)

                !.....transfer from single-lay canopy/soil to soil
                elseif ((pool_indx_lay(n) .eq. 1) .and. &
                        (pool_indx_lay(m) .eq. nsoil)) then
                         do s=1,nsoil
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                              sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                              sibg%l(myl)%poollt%loss_transf_lay(n,1) &
                              * sibg%l(myl)%vegt%rootf(s) &
                              * poolcon(myl)%pool_trans_frac(n,m)
                         enddo
      
                !.....transfer from soil to soil
                elseif ((pool_indx_lay(n) .eq. nsoil) .and. &
                        (pool_indx_lay(m) .eq. nsoil)) then
                         do s=1,nsoil
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                               sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                               sibg%l(myl)%poollt%loss_transf_lay(n,s) &
                               * poolcon(myl)%pool_trans_frac(n,m)
                         enddo
                 else
                      print*, 'Mismatching levels between pool transfers.'
                      print*, 'Stopping in firegfed_interp live transf.'
                      stop
                 endif  !dead pool transfer
              endif !trans_frac > 0.
         enddo  !m=npoolpft/2+1,ntpool/2
     enddo  !n=1,npoolpft/2

     !...same as above but for C13 pools
     do n=npoolpft/2+1,npoolpft !6,10 poolpft
        nref=n+npoolpft/2+1 !12,16 ntpool
        do m=npoolpft/2+1+ntpool/2,ntpool !(17,22) ntpool for C-13 dead pools 
           mref=m-npoolpft !(m-10, 7-12 npoollu dead pools with npoolpft=10)
              if (poolcon(myl)%pool_trans_frac(nref,m) > dzero) then
                !.....transfer from single-layer canopy/soil to surface
                if ((pool_indx_lay(nref) .eq. 1) .and. &
                    (pool_indx_lay(m) .eq. 1)) then
                       sibg%l(myl)%pooldt%gain_transf_lay(mref,1) = &
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,1) + &
                           sibg%l(myl)%poollt%loss_transf_lay(n,1) &
                           * poolcon(myl)%pool_trans_frac(nref,m)

                !.....transfer from single-lay canopy/soil to soil
                elseif ((pool_indx_lay(nref) .eq. 1) .and. &
                        (pool_indx_lay(m) .eq. nsoil)) then
                         do s=1,nsoil
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                             sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                             sibg%l(myl)%poollt%loss_transf_lay(n,1)*sibg%l(myl)%vegt%rootf(s) &
                             * poolcon(myl)%pool_trans_frac(nref,m)
                         enddo

                !.....transfer from soil to soil
                elseif ((pool_indx_lay(nref) .eq. nsoil) .and. &
                        (pool_indx_lay(m) .eq. nsoil)) then
                         do s=1,nsoil
                           sibg%l(myl)%pooldt%gain_transf_lay(mref,s) = &
                             sibg%l(myl)%pooldt%gain_transf_lay(mref,s) + &
                             sibg%l(myl)%poollt%loss_transf_lay(n,s) &
                             * poolcon(myl)%pool_trans_frac(nref,m)
                         enddo
                 else
                      print*, 'Mismatching levels between pool transfers.'
                      print*, 'Stopping in firegfed_interp live transf c13.'
                      stop
                 endif  !dead pool transfer
              endif !trans_frac > 0.
        enddo  !m=npoolpft/2+1+ntpool/2,ntpool
     enddo  !n=npoolpft/2+1,npoolpft


     !-----Carbon Transfers from dead pools-----
     !...j is the sending/from pool
     !...k is the recieving/to pool
     do j=1,npoollu/2 !npoollu=12, (1,6) dead pools
         jref=j+npoolpft/2 ! npoolpft=10, (6,11) from ntpool 
         do k=1,npoollu/2 ! (1,6) dead pools
            kref=k+npoolpft/2 ! (6,11) ntpool
            if (poolcon(myl)%pool_trans_frac(jref,kref) > dzero) then
              !.....transfer gains
              if ((pool_indx_lay(jref) .eq. 1) .and. &
                 (pool_indx_lay(kref) .eq. 1)) then
                      sibg%l(myl)%pooldt%gain_transf_lay(k,1) = &
                        sibg%l(myl)%pooldt%gain_transf_lay(k,1) + &
                        sibg%l(myl)%pooldt%loss_transf_lay(j,1) &
                        * poolcon(myl)%pool_trans_frac(jref,kref)
              elseif ((pool_indx_lay(jref) .eq. 1) .and. &
                     (pool_indx_lay(kref) .eq. nsoil)) then
                      do s=1,nsoil
                        sibg%l(myl)%pooldt%gain_transf_lay(k,s) = &
                           sibg%l(myl)%pooldt%gain_transf_lay(k,s) + &
                           sibg%l(myl)%pooldt%loss_transf_lay(j,1)*sibg%l(myl)%vegt%rootf(s) &
                           * poolcon(myl)%pool_trans_frac(jref,kref)               
                      enddo
              elseif ((pool_indx_lay(jref) .eq. nsoil) .and. &
                     (pool_indx_lay(kref) .eq. nsoil)) then
                      sibg%l(myl)%pooldt%gain_transf_lay(k,:) = &
                        sibg%l(myl)%pooldt%gain_transf_lay(k,:) + &
                        sibg%l(myl)%pooldt%loss_transf_lay(j,:) &
                        * poolcon(myl)%pool_trans_frac(jref,kref)
              else
                   print*,'Mismatching levels between pool transfers.'
                   print*,'Stopping in firegfed_interp dead transf.'
                   stop
              endif !transfer gains
            endif !trans frac > 0
         enddo  !k=1,npoollu
     enddo  !j=1,npoollu

     !...same as above but for c13 pools
     do j=npoollu/2+1,npoollu !(7,12) npoollu dead pools
         jref=j+npoolpft !(17,22) ntpool
         do k=npoollu/2+1,npoollu !7,12 npoollu dead pools
            kref=k+npoolpft !(17,22) ntpool
            if (poolcon(myl)%pool_trans_frac(jref,kref) > dzero) then
             !.....transfer gains
             if ((pool_indx_lay(jref) .eq. 1) .and. &
                 (pool_indx_lay(kref) .eq. 1)) then
                      sibg%l(myl)%pooldt%gain_transf_lay(k,1) = &
                        sibg%l(myl)%pooldt%gain_transf_lay(k,1) + &
                        sibg%l(myl)%pooldt%loss_transf_lay(j,1) &
                        * poolcon(myl)%pool_trans_frac(jref,kref)
             elseif ((pool_indx_lay(jref) .eq. 1) .and. &
                     (pool_indx_lay(kref) .eq. nsoil)) then
                      do s=1,nsoil
                        sibg%l(myl)%pooldt%gain_transf_lay(k,s) = &
                          sibg%l(myl)%pooldt%gain_transf_lay(k,s) + &
                          sibg%l(myl)%pooldt%loss_transf_lay(j,1)*sibg%l(myl)%vegt%rootf(s) &
                          * poolcon(myl)%pool_trans_frac(jref,kref)
                      enddo
             elseif ((pool_indx_lay(jref) .eq. nsoil) .and. &
                     (pool_indx_lay(kref) .eq. nsoil)) then
                      sibg%l(myl)%pooldt%gain_transf_lay(k,:) = &
                        sibg%l(myl)%pooldt%gain_transf_lay(k,:) + &
                        sibg%l(myl)%pooldt%loss_transf_lay(j,:) &
                        * poolcon(myl)%pool_trans_frac(jref,kref)
             else
                   print*,'Mismatching levels between pool transfers.'
                   print*,'Stopping in firegfed_interp dead transf c13.'
                   stop
             endif !transfer gains
            endif !trans frac > 0
          enddo  !k=1,npoollu
      enddo  !j=1,npoollu

      !...save fire respiration
      sibg%l(myl)%poollt%resp_fire = curemis * dtisib
      sibg%l(myl)%poollt%resp_firec13 = curemisc13 * dtisib

     !...Accumulate gains for pool updates
     sibg%l(myl)%pooldt%poollu_dgain = sibg%l(myl)%pooldt%poollu_dgain &!(npoollu,nsoil)
        + sibg%l(myl)%pooldt%gain_transf_lay*dtsib

   enddo !loop thrugh PFTs

!section below is done for individual pools above
!pooldt%poollu_dloss = pooldt%poollu_dloss & !(npoollu,nsoil)
!   + pooldt%loss_resp_lay*dtsib + pooldt%loss_trans_lay*dtsib


      !... Orig. calculation based on GFED4 emissions & rcpoolfire only
      !sibg%l(myl)%poollt%resp_firec13 = tempemisc13 * dtisib
      !... Now updated with actual fire isotope ratio based on pools burned
      !if (tmppooltot .GT. dnzero) then
      !    rcemispoolc13 = (tmppoolc13/tmppooltot)
      !    tempemispoolc13 = rcemispoolc13*tempemis
      !    sibg%l(m yl)%poollt%resp_firec13 = tempemispoolc13 * dtisib
      !    sibg%l(m yl)%fract%d13cemisfire_pool = &
      !       (((1. 0D0/((1.0D0/rcemispoolc13)-1.0D0))/pdb)-1.0D0)*1000.0D0
      !..notes: 1.  post-processing weighting would be done with sum(fire_loss_lay)
      !........ 2.  also in post-processing, assign NaN to this isotope value if
      !........     resp_fire < 0 for any timestep 
      !endif        
                    
   !...Print Result s
   !IF ((fb_err) .O R. (fireb_print)) THEN
   !   print*,''    
   !   print*,'---F IRE CARBON---'
   !   IF (fb_err)  THEN
   !      print('(a )'),'!!Fire Carbon Imbalance!!'
   !      print*,'F ire Emissions Mismatch (mol C/m2): ', curemis
   !   ENDIF        
                    
   !   print('(a,a, i3,a,i4)'), &
   !       '      D ate: ', trim(month_names(month)), day, ', ', year
   !   print('(a,i6 ,2f8.2)'),  '      Point/Lon/Lat: ', indx, lon, lat
   !   print('(a,i14)'),       '      Current Model Second: ', sec_tot
   !   print('(a,2i12)'),      '      Current/Next Fire Second: ', &
   !           fire_seccur, fire_secnext

   !   print*,''
   !   print('(a,i4)'),     '                                ntpft: ', &
   !        ntpft
   !   print*,''
   !   print('(a,f18.8)'),     '                            tfarea: ', &
   !        tfarea

   !   print*,''
   !   print('(a,f18.8)'),     '             sumtotal biomass(m-2): ', &
   !        sum(tpbiomass)
!  !    print*,''
!  !    print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!  !         'biomass(m-2): ', &
!  !         tpbiomass(1),' ',tpbiomass(2),' ',tpbiomass(3),' ',&
!  !         tpbiomass(4),' ',tpbiomass(5)
   !   print*,''
   !   print('(a,f18.8)'),     'sumtotal above-ground biomass(m-2): ', &
   !        sum(tpagb)
!  !    print*,''
!  !    print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!  !         'above-ground biomass(m-2): ', &
!  !         tpagb(1),' ',tpagb(2),' ',tpagb(3),' ',&
!  !         tpagb(4),' ',tpagb(5)

   !   print*,''
   !   print('(a,f18.8)'),     '      Fire C Emissions (umol/m2/s): ', &
   !        sibg%gprogt%firec*mol_to_umol
   !   print('(a,f18.8)'),     '      Time-Step C Losses (mol/m2):', totemis

   !   tempemis = sum(flosspft(:,1:5)) + sum(flosslu(:,1:6))
   !   print('(a,f18.8)'),     '      SiB4 C Removal (mol/m2):    ', tempemis
   !   print('(a)'),           '         PFT    Loss          %-BioBurned   Fraction'
   !      do l=1,ido
   !         myl = tsortref(l)
   !         curemis = sum(flosslu(myl,1:6)) + sum(flosspft(myl,1:5))
   !         if (tpbiomass(myl) .gt. dzero) then
   !            pcemis = curemis/tpbiomass(myl)*100.
   !         else
   !            pcemis = dzero
   !         endif
   !         print('(a,i2,2f14.8,a,f6.2)'), '          ', tpref(myl),  &
   !            curemis, pcemis, '  ',tparea(myl)/tfarea
   !      enddo
   !   print('(a,f12.8)'),     '      Non-Matched C Respired: ', sum(flossb)

   !   IF (fb_err .AND. (fireb_stop)) STOP
   !ENDIF  !print

   deallocate(tpref,tpnum,tparea)
   !deallocate(tpagb,tpagbtemp)
   deallocate(flosspft,flosslu,fgainlu)
   deallocate(flosspft_dm,flosslu_dm)

end subroutine firegfed_interp

