module ic
   use, intrinsic :: iso_c_binding
   use fourier

   implicit none

   namelist /param/ nz, period, lz, lx, nx, nk, nth, delta, a0_peak, xcp, alfa, g_amp, g_x0, g_x1, r0, r1, &
      sigma, gamma, xe, c, x_out, intrvl, in_type, recount, thout, central_mirror, cont, &
      lt, imp_t0, imp_tsp, omega_d, itrec, notk, station, cold, lambda, s2k, norm, ops
    real(c_double)   h, lz, lx, lt, hx, hth, delta, a0_peak, imp_x0, imp_xsp, imp_t0, imp_tsp, g_amp, g_x0, g_x1, r0, r1, x_out, xcp, alfa, &
      sigma, gamma, xe, c, c3, omega_d, kk, lambda, norm
    integer(c_int)   nx, nth, nz, nt, iimp_x0, iimp_xend, iimp_t0, iimp_tend, ifsm_x0, ifsm_end, ix_out, ia_0_x0, ia_0_xend, iasm_x0, iasm_xend, &
      it_out, in_type, nk, intrvl, ig0, ig1, ibe, itrec, ixe1, ixe2
   real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
   logical(c_bool) recount, period, thout, central_mirror, cont, notk, station, cold, s2k, ops

    complex(c_double_complex), allocatable :: a0(:,:), a1(:,:), ak0(:,:), ak1(:,:), jk0(:,:), jk1(:,:), atmp(:), a0z(:), aktmp(:), a0z0cut(:), azl(:), akzl(:), &
                                      dlt(:), ex(:), ak_amp_z0_tend(:), ak_amp_zl_tend(:), akz0(:, :), jkz0(:, :), ak0z0(:), a0z0(:)
   real(c_double), allocatable :: th(:, :, :), dthdz(:, :, :), fk1(:), fk2(:), rhs0(:, :), z(:), x(:), k2(:), &
                                  a_amp_z0_tend(:), a_amp_zl_tend(:), g(:), sum_abs2_a_plus_by_z(:), sum_abs2_a_plus_by_z_k(:)
   complex(c_double_complex), parameter :: im1 = (0.0d0, 1.0d0)
   real(c_double) :: start_time, finish_time

   !interface
   !    subroutine fn_for_fortran_to_call(ptr) &
   !        bind(c, name='fn_for_fortran_to_call')
   !        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
   !        implicit none
   !        type(c_ptr), intent(in), value :: ptr
   !        !integer(c_int), intent(in), value :: arg
   !    end subroutine fn_for_fortran_to_call
   !end interface

contains
   subroutine calc_idx()
      implicit none

      if (lt .lt. 2) then
         print *, 'lt must be greater or equal then 2'
         pause
         stop
      end if

      kk = 4.0d0*pi/lambda
      if (s2k .eq. .true.) then
         !lx = lx / sqrt(kk)

         !what depends on the new lx
         xcp = 0.5d0*lx! * norm ** (1.0/6.0)
         alfa = 0.5d0*lx/100.0d0! * norm ** (1.0/6.0)
         g_x0 = 0.35*lx! * norm ** (1.0/6.0)
         g_x1 = 0.65*lx! * norm ** (1.0/6.0)
         !xe = 0.166666666666667 * lx! * norm ** (1.0/6.0)
         xe = 0.22*lx! * norm ** (1.0/6.0)
         x_out = 0.5d0*lx! * norm ** (1.0/6.0)

         s2k = .false.
      end if

      if (norm .ne. 0.0d0) then
         lx = lx*norm**(1.0/6.0)

         !what depends on the new lx
         xcp = 0.5d0*lx! * norm ** (1.0/6.0)
         alfa = 0.5d0*lx/100.0d0! * norm ** (1.0/6.0)
         g_x0 = 0.35*lx! * norm ** (1.0/6.0)
         g_x1 = 0.65*lx! * norm ** (1.0/6.0)
         !xe = 0.166666666666667 * lx! * norm ** (1.0/6.0)
         xe = 0.22*lx! * norm ** (1.0/6.0)
         x_out = 0.5d0*lx! * norm ** (1.0/6.0)

         delta = delta/norm**(1.0/3.0)
         a0_peak = a0_peak/norm**(2.0/3.0)
         sigma = sigma/norm**(1.0/3.0)
         c = c/norm
         if (period == .true.) then
            lz = lz*(lx*lx)/lambda
         end if
         period = .false.
         norm = 0
      else
         if (period == .true.) then
            lz = lz*(lx*lx)/lambda
         end if
         period = .false.
      end if

      open (unit=1, file='input_fortran_real.in')
      write (unit=1, nml=param)
      close (unit=1)

      !c3 = c * c * c
      c3 = c

      h = lz/nz
      nt = nz*lt + 1 !lt - in shares of lz
      nz = nz + 1

      hth = 2.0d0*pi/nth; 
      hx = lx/nx; 
      ixe1 = xe/hx + 1
      ixe2 = nx - ixe1 + 2
      xe = (ixe1 - 1)*hx !xe update

      print *, 'ixe1 = ', ixe1
      print *, 'ixe2 = ', ixe2
      print *, ' xe = ', xe
      print *, '(ixe1 - 1) * hx = ', (ixe1 - 1)*hx
      print *, ' lx - xe = ', lx - xe
      print *, '(ixe2 - 1) * hx = ', (ixe2 - 1)*hx
      !if ((lx-xe) == (ixe2 - 1) * hx) print *, '(lx-xe) == (ixe2 - 1) * hx'

      !beginning of the end
      ibe = nt - nz + 1
      if (ibe .le. 0) then
         print *, 'ibe must be greater 0'
         pause
         stop
      end if

      ix_out = int(x_out/hx)
      if (ix_out <= 0) ix_out = 1
      if (ix_out > nx) ix_out = nx

      it_out = intrvl

      iimp_x0 = max(1, int(imp_x0/hx) + 1)
      imp_x0 = (iimp_x0 - 1)*hx !to (iimp_x0 - 1) * hx exactly equal to imp_x0
      iimp_xend = min(nx + 1, int((imp_x0 + imp_xsp)/hx) + 1) ! calculate for nx'= nx + 1
      imp_xsp = (iimp_xend - iimp_x0)*hx !for hit precision (-1 + 1 = 0)
      if (iimp_xend == nx + 1) iimp_xend = iimp_xend - 1 ! last interval point not involved

      iimp_t0 = max(1, int(imp_t0/h) + 1)
      iimp_tend = nz*imp_tsp + iimp_t0 !imp_tsp - in shares of lz
      !imp_t0 = (iimp_t0 - 1) * h; !to (iimp_t0 - 1) * h exactly equals imp_t0
      !iimp_tend = min(nt, int((imp_t0 + imp_tsp) / h) + 1) ! counting for nt
      imp_tsp = (iimp_tend - iimp_t0)*h !for point hit accuracy (-1 + 1 = 0) (imp_tsp - no longer in shares of lz) !!!

      iasm_x0 = max(1, int(g_x0/hx) + 1)
      g_x0 = (iasm_x0 - 1)*hx
      iasm_xend = min(nx + 1, int((g_x0 + g_x1)/hx) + 1) ! calculate for nx'= nx + 1
      g_x1 = (iasm_xend - iasm_x0)*hx
      if (iasm_xend == nx + 1) iasm_xend = iasm_xend - 1 ! last interval point not involved
   end subroutine

   subroutine calc_theta(th, dthdz)
      implicit none

      real(c_double), intent(inout) :: th(:, :), dthdz(:, :)
      integer, dimension(size(th)) :: i
      integer ix

      i = (/1:size(th, 1)/)

      do ix = 1, 2
         th(:, ix) = hth*(i - 1)
         dthdz(:, ix) = delta
      end do
   end subroutine
end module ic

program elektron2d
   !subroutine calculate_fortran(ptr) bind(c,name='calculate_fortran')
   use, intrinsic :: iso_c_binding, only: c_ptr
   use ic, only: h, hx, th, dthdz, rhs0, a0, a1, nz, nx, nt, nth, recount, atmp, r0, r1, g, im1, k2, lz, &
                 ak0, ak1, jk0, jk1, delta, ix_out, a_amp_z0_tend, a_amp_zl_tend, &
                 ak_amp_z0_tend, ak_amp_zl_tend, start_time, finish_time, calc_theta, z, x, &
                 it_out, intrvl, fk1, fk2, nk, dlt, ex, gamma, sigma, c3, &
 ig0, ig1, sum_abs2_a_plus_by_z, sum_abs2_a_plus_by_z_k, lx, cont, omega_d, akz0, jkz0, itrec, notk, ixe1, ixe2, xe, station, a0z, &
                 ak0z0, a0z0, akzl, a0z0cut, azl, aktmp, c
   use fourier
   use ifport

   implicit none

   interface
      subroutine init() bind(c, name='init')
      end subroutine init
      subroutine finish() bind(c, name='finish')
      end subroutine finish
      function a0_fn(it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nx) :: a
         integer(c_int), intent(in) :: it
      end function a0_fn
      function a0_fn_stat(n_it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nx) :: a
         integer(c_int), intent(in) :: n_it
      end function a0_fn_stat
      function ak0_fn(it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nk) :: a
         integer(c_int), intent(in) :: it
      end function ak0_fn
      function ak0_fn_stat(n_it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nk) :: a
         integer(c_int), intent(in) :: n_it
      end function ak0_fn_stat
      subroutine makea(atmp, ak)
         use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
         complex(c_double_complex), dimension(:), intent(inout) :: atmp
         complex(c_double_complex), dimension(:), intent(in) :: ak
      end subroutine makea
      subroutine makeak(ak, atmp)
         use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
         complex(c_double_complex), dimension(:), intent(inout) :: ak
         complex(c_double_complex), dimension(:), intent(in) :: atmp
      end subroutine makeak
      subroutine make_a0z0_ak1_atmp(a0z0, ak1, atmp, ak)
         use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
         complex(c_double_complex), dimension(:), intent(inout) :: a0z0, ak1, atmp
         complex(c_double_complex), dimension(:), intent(in) :: ak
      end subroutine make_a0z0_ak1_atmp
   end interface

   !type(c_ptr), intent(in), value :: ptr
   integer(c_int) :: iz, ix, it, percent, i, j, k, itau0, rec_length, nr, nzap = 0
   integer(c_int) :: istat, int_time, ires
    real(c_double) float_time, eff_tmp, eff_tmp_k, sum_abs2_a_by_zx_plus, sum_abs2_ak_by_zx_plus, sum_abs2_a_by_zx_minus, sum_abs2_ak_by_zx_minus, &
        int_abs2_a_plus, int_abs2_a_minus, eff(2), int_abs2_a_minus_at_zl_on_mir, int_abs2_a_minus_at_zl_on_mir_k, loss_on_the_way_minus_k, int_abs2_a_minus_at_z0_on_mir, &
        int_abs2_a_minus_at_z0_out_mir, loss_on_the_way_minus, int_abs2_a_plus_at_z0, int_abs2_a_plus_at_z0_k, omega, sum_eff, int_abs2_a_plus_at_zl, &
      int_abs2_a_plus_at_zl_k, int_abs2_a_plus_at_zl_on_mir, int_abs2_a_plus_at_zl_on_mir_k, &
      int_abs2_a_plus_at_zl_out_mir, int_abs2_a_plus_at_zl_out_mir_k, loss_on_the_way_plus, loss_on_the_way_plus_k, omega_dd, &
        int_abs2_a_minus_at_z0_on_mir_k, int_abs2_a_minus_at_z0_out_mir_k, int_abs2_a_minus_at_z0_k, int_abs2_a_minus_at_z0, int_obrezki_z0, &
      int_obrezki_zl, obrezki, eff_tmp_b, eff_tmp_k_b
   complex(c_double_complex) :: aend0 = 0, aend1 = 0
   character*10 str

   call init()

   !what does not depend on j
   dlt = im1*gamma*k2 + sigma
   ex = cdexp(-dlt*h)

   !at the beginning:
   jk0 = 0; th = 0; dthdz = 0 ! for (z=0,tau=0)
   call calc_theta(th(:, 1, :), dthdz(:, 1, :)) ! for z=0 for all tau (for the other z, with tau=0, not used)
   jk1(1, :) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2)) !initial conditions for current at z=0

   !for zero layer:
   if (recount == .false.) then
      call makea(atmp, ak0_fn_stat(0))

      atmp(:) = atmp*r0*g !and cut

      call makeak(ak0(1, :), atmp)
      jk0(1, :) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2))
   else
      if (cont == .true.) then
         open (2, file='akzxt0.bin', form='binary', err=101)
         open (3, file='jkzxt0.bin', form='binary', err=101)
         do iz = 1, nz
            read (2, err=102) ak0(iz, :)
            read (3, err=102) jk0(iz, :)
         end do
         close (2)
         close (3)
         do iz = 1, nz
            ak0(iz, :) = ak0(iz, :)*cdexp(-im1*omega_d*(iz - 1)*h)
            jk0(iz, :) = jk0(iz, :)*cdexp(-im1*omega_d*(iz - 1)*h)
         end do
      else
         call akt0_fn(ak0) !initial conditions for tau=0
         do i = 1, nz
            jk0(i, :) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2))
         end do
         open (1, file='akzxt0.bin', form='binary', err=101)
         open (3, file='jkzxt0.bin', form='binary', err=101)
         do iz = 1, nz
            write (1, err=103) ak0(iz, :)
            write (3, err=103) jk0(iz, :)
         end do
         close (1)
         close (3)
         open (1, file='ak_t0.bin', form='binary', err=101)
         open (3, file='jk_t0.bin', form='binary', err=101)
         do iz = 1, nz
            write (1, err=103) ak0(iz, :)
            write (3, err=103) jk0(iz, :)
         end do
         close (1)
         close (3)
         open (1, file='azxt0.bin', form='binary', err=101)
         do iz = 1, nz
            call makea(atmp(:), ak0(iz, :))
            write (1, err=103) atmp
         end do
         close (1)
      end if
   end if

!for z=0:
   if (recount == .false.) then
      do it = 1, nz
         akz0(:, it) = ak0_fn_stat(0)
         jkz0(:, it) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2))
      end do
   else
      if (cont == .true.) then
         if (station .eq. .true.) then
            do it = 1, nz
               akz0(:, it) = ak0(1, :)*cdexp(im1*omega_d*(it - 1)*h)
               jkz0(:, it) = jk0(1, :)*cdexp(im1*omega_d*(it - 1)*h)
            end do
         else
            open (1, file='akxtz0.bin', form='binary', err=101)
            open (3, file='jkxtz0.bin', form='binary', err=101)
            do it = 1, nz
               read (1, err=102) akz0(:, it)
               read (3, err=102) jkz0(:, it)
            end do
            close (1)
            close (3)
            do it = 1, nz
               akz0(:, it) = akz0(:, it)*cdexp(im1*omega_d*(it - 1)*h)
               jkz0(:, it) = jkz0(:, it)*cdexp(im1*omega_d*(it - 1)*h)
            end do
         end if
      else
         do it = 1, nz
            akz0(:, it) = ak0_fn(it)
            jkz0(:, it) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2))
         end do
         open (355, file='ak_t.bin', form='binary', buffered='no', err=101)
         open (356, file='jk_t.bin', form='binary', buffered='no', err=101)
         do it = 1, nz
            write (355) akz0(:, it)
            write (356) jkz0(:, it)
         end do
         close (355)
         close (356)
      end if
   end if

   !remember for history
   ak0z0(:) = akz0(:, 1)

   !ak0(1,:) == akz0(:,1)

   !and cut it
   call makea(atmp(:), ak0(1, :))
   atmp(:) = atmp*g !cutting
   call makeak(ak0(1, :), atmp) !and initialize

   !record length
   rec_length = c_double_complex*2*nk/4
   nr = 1
   !open files
   open (222, file='$$$.bin', access='direct', recl=rec_length, status='replace', buffered='no', err=101)
   open (8, file='t.dat', buffered='no', err=101)
   open (3, file='eff.dat', err=101)
   open (35, file='eff_b.dat', err=101)
   open (53, file='eff_new.dat', err=101)
   open (33, file='eff_k.dat', err=101)
   open (335, file='eff_k_b.dat', err=101)
   open (533, file='eff_new_k.dat', err=101)
   open (553, file='ari.dat', buffered='no', err=101)
   open (777, file='for_graphics.dat', err=101)
   !we write what is at the time tau=0
   write (222, rec=nr) ak0(nz, :)
   !print *, ' initial nr = ', nr

   write (*, '(/)')

   do_t: do it = 1, nt
      if (itrec .eq. 0 .or. mod(it, (nz - 1)*itrec) .eq. 1) then
         open (1, file='nit_tmp.dat', err=101)
         write (1, *) 'at the beginning of the last block:'
         write (1, *) 'it = ', it
         write (1, *) 'tau = ', (it - 1)*h
         write (1, *) 'omega = ', omega
         write (1, *) 'delta = ', delta
         write (1, *) 'ixe1 = ', ixe1
         write (1, *) 'ixe2 = ', ixe2
         write (1, *) 'xe = ', xe
         write (1, *) '(ixe1 - 1) * hx = ', (ixe1 - 1)*hx
         write (1, *) '(ixe2 - 1) * hx = ', (ixe2 - 1)*hx
         write (1, *) 'ht = h = ', h
         close (1)
         if (nzap .eq. 0) nzap = 1 ! countdown
         open (880, file='akxtz0_tmp.bin', form='binary', buffered='no', err=101)
         open (881, file='akxtzl_tmp.bin', form='binary', buffered='no', err=101)
         open (882, file='jkxtz0_tmp.bin', form='binary', buffered='no', err=101)
      end if

      if (nzap .eq. nz) then
         open (888, file='akxtz0_n_tmp.bin', form='binary', buffered='no', err=101)
         open (889, file='akxtzl_n_tmp.bin', form='binary', buffered='no', err=101)
         open (890, file='jkxtz0_n_tmp.bin', form='binary', buffered='no', err=101)
      end if

      if (nzap .ge. 1 .and. nzap .le. nz) then
         !a0 write before cutting
         write (880) ak0z0
         write (881) ak0(nz, :)
         write (882) jk0(1, :)
      end if

      if (nzap .ge. nz .and. nzap .le. 2*nz - 1) then
         !a0 write before cutting
         write (888) ak0z0
         write (889) ak0(nz, :)
         write (890) jk0(1, :)
      end if

      if (nzap .eq. 1) then
         open (883, file='azxt0_tmp.bin', form='binary', buffered='no', err=101)
         open (886, file='akzxt0_tmp.bin', form='binary', buffered='no', err=101)
         open (884, file='jkzxt0_tmp.bin', form='binary', buffered='no', err=101)
         call makea(a0z0, ak0z0)
         write (883, err=103) a0z0
         write (886, err=103) ak0z0
         write (884, err=103) jk0(1, :)
      end if

      if (nzap .eq. nz) then
         open (885, file='azxtl_tmp.bin', form='binary', buffered='no', err=101)
         open (887, file='akzxtl_tmp.bin', form='binary', buffered='no', err=101)
         call makea(a0z0, ak0z0)
         write (885, err=103) a0z0
         write (887, err=103) ak0z0
      end if

      if (nzap .eq. nz + 1) then
         open (892, file='akzxt0_m_tmp.bin', form='binary', buffered='no', err=101)
         open (893, file='jkzxt0_m_tmp.bin', form='binary', buffered='no', err=101)
         write (892, err=103) ak0z0
         write (893, err=103) jk0(1, :)
      end if

      if (nzap .eq. 2*nz - 1) then
         open (895, file='akzxtl_n_tmp.bin', form='binary', buffered='no', err=101)
         write (895, err=103) ak0z0
      end if

      !initial conditions for theta and detune at z=0
      call calc_theta(th(:, 1, :), dthdz(:, 1, :))

      !initial conditions for current at z=0
      jk1(1, :) = fk1*jf(th(:, 1, 1)) + fk2*jf(th(:, 1, 2))

      if (recount .eq. .true.) then
         if (it >= nz - 1) then
            itau0 = nr - nz + 2

            akzl = dcmplx(0)
            read (222, rec=itau0) akzl !all what comes in zl is written in k-space

            call makea(atmp, akzl) !make uncut field on lz

            atmp = r1*atmp*g !mirror reflection and cut to z=lz

            !fourier
            call makeak(akzl, atmp) !cut reflected field and cut harmonics from it
            call makea(atmp, akzl) !cut harmonic field

            int_abs2_a_minus_at_zl_on_mir_k = 0.5d0*dmysum(cdabs(akzl)*cdabs(akzl)) !for efficiency (that is already reflected from the mirror) (k-space)
            int_abs2_a_minus_at_zl_on_mir = sum(cdabs(atmp)*cdabs(atmp))*hx ! in x-space

            !-------------------------------------------------------------------------------------------------------------------------------------------------------

            ak0z0 = akzl*cdexp(-dlt*lz) !that is back to z0, (k-space)

            call make_a0z0_ak1_atmp(a0z0, a0z0cut, ak1(1, :), ak0z0)

            int_abs2_a_minus_at_z0_k = 0.5d0*dmysum(cdabs(ak0z0)*cdabs(ak0z0)) !energy of what is back in z0 (k-space)
            int_abs2_a_minus_at_z0 = sum(cdabs(a0z0)*cdabs(a0z0))*hx !energy of what is back in z0 (x-space)

            int_abs2_a_minus_at_z0_on_mir = sum(cdabs(a0z0(ig0:ig1 - 1))*cdabs(a0z0(ig0:ig1 - 1)))*hx !(x-space)
                int_abs2_a_minus_at_z0_out_mir = (sum(cdabs(a0z0(1:ig0-1))*cdabs(a0z0(1:ig0-1))) + sum(cdabs(a0z0(ig1+1:nx))*cdabs(a0z0(ig1+ 1:nx)))) * hx !(x-space)

            loss_on_the_way_minus_k = int_abs2_a_minus_at_zl_on_mir_k - int_abs2_a_minus_at_z0_k
            loss_on_the_way_minus = int_abs2_a_minus_at_zl_on_mir - (int_abs2_a_minus_at_z0_on_mir + int_abs2_a_minus_at_z0_out_mir) ! subtract what is returned to the beginning, the rest - loss (x-space)

            aktmp = a0z0*g
            call sint(aktmp) ! uncut reflection cut modes
            int_abs2_a_minus_at_z0_out_mir_k = int_abs2_a_minus_at_z0_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

            aktmp(1:nk) = dcmplx(0) ! cut mode of reflection cut
            int_obrezki_z0 = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to z0

                !!fourier
            !call makeak(akzl, atmp)
            !
            !int_abs2_a_minus_at_zl_on_mir_k = 0.5d0 * dmysum(cdabs(akzl) * cdabs(akzl))
            !
            !call make_a0z0_ak1_atmp(a0z0, ak1(1,:), a0z0cut, ak0z0)
            !
            !open(388, file = 'a0.bin', form = 'binary', err = 101)
            !write(388) a0z0cut
            !close(388)
            !
            !to_chto_vernulos_v_nachalo_k = 0.5d0 * dmysum(cdabs(ak0z0) * cdabs(ak0z0))
            !loss_on_the_way_minus_k = int_abs2_a_minus_at_zl_on_mir_k - &
            ! to_chto_vernulos_v_nachalo_k
            !
                !! the same losses as before only in x-space
            !int_abs2_a_minus_at_z0_on_mir = sum(cdabs(a0z0(ig0:ig1-1))*cdabs(a0z0(ig0:ig1-1))) * hx
            !int_abs2_a_minus_at_z0_out_mir = (sum(cdabs(a0z0(1:ig0-1))*cdabs(a0z0(1:ig0-1))) + sum(cdabs(a0z0(ig1+1:nx))*cdabs(a0z0(ig1 +1:nx)))) * hx
            !loss_on_the_way_minus = int_abs2_a_minus_at_zl_on_mir - (int_abs2_a_minus_at_z0_on_mir + int_abs2_a_minus_at_z0_out_mir)
            !
            !int_abs2_a_minus_at_z0_on_mir_k = 0.5d0 * dmysum(cdabs(ak1(1,:)) * cdabs(ak1(1,:)))
            !int_abs2_a_minus_at_z0_out_mir_k = to_chto_vernulos_v_nachalo_k - int_abs2_a_minus_at_z0_on_mir_k
         else
            ak0z0 = akz0(:, it + 1)

            call make_a0z0_ak1_atmp(a0z0, a0z0cut, ak1(1, :), ak0z0)
         end if
      else
         ak0z0 = ak0_fn_stat(0)

         call make_a0z0_ak1_atmp(a0z0, a0z0cut, ak1(1, :), ak0z0)
      end if

      !to calculate the efficiency we calculate the sum of abs(a(z=0))**2 by x
      int_abs2_a_plus_at_z0 = sum(cdabs(a0z0cut)*cdabs(a0z0cut))*hx !(x-space)
      int_abs2_a_plus_at_z0_k = 0.5d0*sum(cdabs(ak1(1, :))*cdabs(ak1(1, :))) !(k-space)

      !initial sum of square field amplitudes at z=0 (also for efficiency) (vector)
      sum_abs2_a_plus_by_z = cdabs(a0z0cut)*cdabs(a0z0cut) !(x-space)
      sum_abs2_a_plus_by_z_k = 0.5d0*cdabs(ak1(1, :))*cdabs(ak1(1, :)) !(k-space)

      do_z: do iz = 1, nz - 1
         rhs0 = rhs(ak1(iz, :), th(:, iz, :))
         th(:, iz + 1, :) = th(:, iz, :) + dthdz(:, iz, :)*h + h/2.0d0*rhs0*h !theta predictor
         jk1(iz + 1, :) = fk1*jf(th(:, iz + 1, 1)) + fk2*jf(th(:, iz + 1, 2)) !current predictor

         !predictor a (interpolation)
         ak1(iz + 1, 1) = ak0(iz, 1) + h/2.0d0*(jk0(iz, 1) + jk1(iz + 1, 1))
         ak1(iz + 1, 2:nk) = ak0(iz, 2:nk)*ex(2:nk) + &
                             c3*(jk0(iz, 2:nk) + jk1(iz + 1, 2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                                 ex(2:nk)*(jk1(iz + 1, 2:nk) - jk0(iz, 2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

         !predictor a (keystone)
         !atmp = (ak0(iz,:) + c3 * h / 2.0d0 * jk0(iz,:)) * cdexp(-dlt * h) !part a
         !ak1(iz+1,:) = atmp + c3 * h / 2.0d0 * jk1(iz+1,:) !predictor a

         !theta corrector
         th(:, iz + 1, :) = th(:, iz, :) + dthdz(:, iz, :)*h + h/6.0d0*rhs0*h &
                            + h/3.0d0*rhs(ak1(iz + 1, :), th(:, iz + 1, :))*h

         !j1(iz+1,:) = j_fn(th(:,iz+1,:)) !current corrector
         jk1(iz + 1, :) = fk1*jf(th(:, iz + 1, 1)) + fk2*jf(th(:, iz + 1, 2)) !current corrector

         !offset a (interpolation)
         ak1(iz + 1, 1) = ak0(iz, 1) + h/2.0d0*(jk0(iz, 1) + jk1(iz + 1, 1))
         ak1(iz + 1, 2:nk) = ak0(iz, 2:nk)*ex(2:nk) + &
                             c3*(jk0(iz, 2:nk) + jk1(iz + 1, 2:nk)*(-1.0d0 + dlt(2:nk)*h) + &
                                 ex(2:nk)*(jk1(iz + 1, 2:nk) - jk0(iz, 2:nk)*(1.0d0 + dlt(2:nk)*h)))/dlt(2:nk)/dlt(2:nk)/h

         !corrector a (keystone)
         !atmp = (ak0(iz,:) + c3 * h / 2.0d0 * jk0(iz,:)) * cdexp(-dlt * h) !part a
         !ak1(iz+1,:) = atmp + c3 * h / 2.0d0 * jk1(iz+1,:)

         dthdz(:, iz + 1, :) = dthdz(:, iz, :) + h/2.0d0*(rhs0 + rhs(ak1(iz + 1, :), th(:, iz + 1, :)))

         !calculation of efficiency at point z = iz*h
         eff(1) = 1.0d0/nth*sum(dthdz(:, iz + 1, 1) - delta) !counting efficiency (xcp1)
         eff(2) = 1.0d0/nth*sum(dthdz(:, iz + 1, 2) - delta) !counting efficiency (xcp2)

         !back to reality
         call makea(atmp(:), ak1(iz + 1, :))

         aend1 = atmp(ix_out) !to calculate omega

         !sum of field amplitudes at z=iz*h
         sum_abs2_a_plus_by_z = sum_abs2_a_plus_by_z + cdabs(atmp)*cdabs(atmp)
         sum_abs2_a_plus_by_z_k = sum_abs2_a_plus_by_z_k + 0.5d0*cdabs(ak1(iz + 1, :))*cdabs(ak1(iz + 1, :))

         if (nzap .eq. 1) then
            call makea(a0z, ak0(iz + 1, :))
            write (883, err=103) a0z
            write (886, err=103) ak0(iz + 1, :)
            write (884, err=103) jk0(iz + 1, :)
         end if

         if (nzap .eq. nz) then
            call makea(a0z, ak0(iz + 1, :))
            write (885, err=103) a0z
            write (887, err=103) ak0(iz + 1, :)
         end if

         if (nzap .eq. nz + 1) then
            write (892, err=103) ak0(iz + 1, :)
            write (893, err=103) jk0(iz + 1, :)
         end if

         if (nzap .eq. 2*nz - 1) then
            write (895, err=103) ak0(iz + 1, :)
         end if
      end do do_z

      !open(21, file='test.dat')
      !do i=1,nx
      !    write(21, '(2e17.8)') (i-1)*hx, cdabs(atmp(i))
      !enddo
      !close(21)
      !stop

      if (nzap .ne. 0) nzap = nzap + 1
      if (nzap .gt. 2*nz - 1) then
         close (1)
         close (883) !azxt0.bin
         close (885) !azxtl.bin
         close (880) !akxtz0_tmp.bin
         close (881) !akxtzl_tmp.bin
         close (882) !jkxtz0_tmp.bin
         close (884) !jkzxt0_tmp.bin
         close (886) !akzxt0_tmp.bin
         close (887) !akzxtl_tmp.bin
         close (888) !axtz0_n_tmp.bin
         close (889) !axtzl_n_tmp.bin
         close (890) !jkxtz0_n_tmp.bin
         close (892) !akzxt0_m_tmp.bin
         close (893) !jkzxt0_m_tmp.bin
         close (895) !akzxtl_n_tmp.bin

         nzap = 0

         ires = delfilesqq('nit.dat')
         ires = delfilesqq('axtz0.bin')
         ires = delfilesqq('axtzl.bin')
         ires = delfilesqq('jkxtz0.bin')

         ires = renamefileqq('nit_tmp.dat', 'nit.dat')
         ires = renamefileqq('axtz0_tmp.bin', 'axtz0.bin')
         ires = renamefileqq('axtzl_tmp.bin', 'axtzl.bin')
         ires = renamefileqq('jkxtz0_tmp.bin', 'jkxtz0.bin')

         ires = delfilesqq('axtz0_n.bin')
         ires = delfilesqq('axtzl_n.bin')
         ires = delfilesqq('jkxtz0_n.bin')

         ires = renamefileqq('axtz0_n_tmp.bin', 'axtz0_n.bin')
         ires = renamefileqq('axtzl_n_tmp.bin', 'axtzl_n.bin')
         ires = renamefileqq('jkxtz0_n_tmp.bin', 'jkxtz0_n.bin')

         ires = delfilesqq('azxt0.bin')
         ires = delfilesqq('jkzxt0.bin')
         ires = delfilesqq('azxtl.bin')
         ires = delfilesqq('akzxt0.bin')
         ires = delfilesqq('akzxtl.bin')

         ires = renamefileqq('azxt0_tmp.bin', 'azxt0.bin')
         ires = renamefileqq('jkzxt0_tmp.bin', 'jkzxt0.bin')
         ires = renamefileqq('azxtl_tmp.bin', 'azxtl.bin')
         ires = renamefileqq('akzxt0_tmp.bin', 'akzxt0.bin')
         ires = renamefileqq('akzxtl_tmp.bin', 'akzxtl.bin')

         ires = delfilesqq('azxt0_m.bin')
         ires = delfilesqq('jkzxt0_m.bin')
         ires = delfilesqq('akzxt0_m.bin')

         ires = renamefileqq('azxt0_m_tmp.bin', 'azxt0_m.bin')
         ires = renamefileqq('jkzxt0_m_tmp.bin', 'jkzxt0_m.bin')
         ires = renamefileqq('akzxt0_m_tmp.bin', 'akzxt0_m.bin')

         ires = delfilesqq('azxtl_n.bin')
         ires = delfilesqq('akzxtl_n.bin')

         ires = renamefileqq('azxtl_n_tmp.bin', 'azxtl_n.bin')
         ires = renamefileqq('akzxtl_n_tmp.bin', 'akzxtl_n.bin')
      end if

      !phase
      !omega = dimag(cdlog(a1(nz,ix_out) / a0(nz,ix_out))) / h
      omega = dimag(cdlog(aend1/aend0))/h
      aend0 = aend1 !swap

      !calculation of the sum of efficiency for x at the point z=lz
      sum_eff = (eff(1) + eff(2))/2.0d0
      !to calculate the efficiency we calculate the sum of abs(a(z=lz))**2 by x
      int_abs2_a_plus_at_zl = sum(cdabs(atmp)*cdabs(atmp))*hx
      int_abs2_a_plus_at_zl_k = 0.5d0*dmysum(cdabs(ak1(nz, :))*cdabs(ak1(nz, :)))

      int_abs2_a_plus_at_zl_on_mir = sum(cdabs(atmp(ig0:ig1 - 1))*cdabs(atmp(ig0:ig1 - 1)))*hx
        int_abs2_a_plus_at_zl_out_mir = (sum(cdabs(atmp(1:ig0-1))*cdabs(atmp(1:ig0-1))) + sum(cdabs(atmp(ig1+1:nx))*cdabs(atmp(ig1+1:nx)))) * hx

      aktmp(:) = atmp*g
      call sint(aktmp)
      int_abs2_a_plus_at_zl_out_mir_k = int_abs2_a_plus_at_zl_k - 0.5d0*sum(cdabs(aktmp)*cdabs(aktmp))

      aktmp(1:nk) = dcmplx(0)
      int_obrezki_zl = 0.5d0*dmysum(cdabs(aktmp)*cdabs(aktmp)) ! energy in cut to zl

      !double sum of field amplitudes in z and in x
      loss_on_the_way_plus = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z)*hx*h !x-space
      loss_on_the_way_plus_k = 2.0d0*sigma*sum(sum_abs2_a_plus_by_z_k)*h !k-space

      !azl = r1 * atmp * g !mirror reflection in z=lz
      !call makeak(akzl, azl) !make harmonics from cuttings
      !int_abs2_a_plus_at_zl_on_mir_k = 0.5d0 * dmysum(cdabs(akzl) * cdabs(akzl))
      !
      !int_abs2_a_plus_at_zl_out_mir_k = int_abs2_a_plus_at_zl_k - int_abs2_a_plus_at_zl_on_mir_k

      !recording efficiency at this iteration
      eff_tmp = (loss_on_the_way_plus + int_abs2_a_plus_at_zl - int_abs2_a_plus_at_z0)/lx - 4.0d0*c3*sum_eff
      eff_tmp_k = loss_on_the_way_plus_k + int_abs2_a_plus_at_zl_k - int_abs2_a_plus_at_z0_k - 4.0d0*c3*sum_eff

      obrezki = int_obrezki_z0 + int_obrezki_zl

      write (553, '(4e17.8)', err=103) it*h, dreal(aend1), dimag(aend1), cdabs(aend1)

      if (((intrvl > 0) .and. (mod(it, it_out) == 0)) .or. (it == nt - 1)) then
         write (3, 104, err=103) & ! eff.dat
            it, &
            int_abs2_a_plus_at_zl/lx, &
            loss_on_the_way_plus/lx, &
            int_abs2_a_plus_at_z0/lx, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp, &
            eff_tmp/(4.0d0*c3*sum_eff), &
            omega

         eff_tmp_b = (loss_on_the_way_minus + int_abs2_a_minus_at_z0 - int_abs2_a_minus_at_zl_on_mir)/lx

         write (35, 108, err=103) & ! eff_b.dat
            it, &
            loss_on_the_way_minus/lx, &
            int_abs2_a_minus_at_zl_on_mir/lx, &
            int_abs2_a_minus_at_z0/lx, &
            eff_tmp_b

108      format(i, 4f17.8)

         write (33, 104, err=103) & ! eff_k.dat
            it, &
            int_abs2_a_plus_at_zl_k, &
            loss_on_the_way_plus_k, &
            int_abs2_a_plus_at_z0_k, &
            sum_eff, &
            int_obrezki_z0, &
            int_obrezki_zl, &
            eff_tmp_k, &
            eff_tmp_k/(4.0d0*c3*sum_eff), &
            omega

104      format(i, 9f17.8)

         eff_tmp_k_b = loss_on_the_way_minus_k + int_abs2_a_minus_at_z0_k - int_abs2_a_minus_at_zl_on_mir_k

         write (335, 108, err=103) & ! eff_k_b.dat
            it, &
            loss_on_the_way_minus_k, &
            int_abs2_a_minus_at_zl_on_mir_k, &
            int_abs2_a_minus_at_z0_k, &
            eff_tmp_k_b

         write (53, 107, err=103) it, &   ! eff_new.dat
            int_abs2_a_minus_at_z0_out_mir/lx, &
            int_abs2_a_plus_at_zl_out_mir/lx, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus/lx), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                       loss_on_the_way_plus + &
                                       int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir + &
                                        loss_on_the_way_plus + &
                                        int_abs2_a_minus_at_z0_out_mir)/lx + obrezki - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

         write (533, 107, err=103) it, &   ! eff_new_k.dat
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                       loss_on_the_way_plus_k + &
                                       int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
            4.0d0*c3*sum_eff, &
            (loss_on_the_way_minus_k + (int_abs2_a_plus_at_zl_out_mir_k + &
                                        loss_on_the_way_plus_k + &
                                        int_abs2_a_minus_at_z0_out_mir_k + obrezki) - &
             4.0d0*c3*sum_eff)/4.0d0/c3/sum_eff

107      format(i, 7f17.8)

         write (777, '(9e17.8)', err=103) & ! for_graphics.dat
            lz, &
            int_abs2_a_minus_at_z0_out_mir_k, &
            int_abs2_a_plus_at_zl_out_mir_k, &
            (loss_on_the_way_minus_k + loss_on_the_way_plus_k), &
            obrezki, &
            4.0d0*c3*sum_eff, &
            (int_abs2_a_minus_at_z0_out_mir_k + int_abs2_a_plus_at_zl_out_mir_k)/(4.0d0*c3*sum_eff), &
            omega, &
            sum_eff/c
      end if

        !!dir$ if defined (testflag)
      !if ((recount .eq. .false.) .or. (station .eq. .true.)) then
      if ((recount .eq. .false.)) then
         !a1(1,:) = a0(1,:)
         !j1(1,:) = j0(1,:)
         if (it <= nz) then
            do i = 1, it
               do j = 1, nk
                  if (ak1(i, j) /= ak0(i, j)) then
                     write (*, '(/)')
                     print *, 'z = ', i
                     print *, 'x = ', j
                     print *, 't = ', it
                     write (*, '(a,i3,i3,a,2e17.8)') 'ak1(', i, j, ') = ', ak1(i, j)
                     write (*, '(a,i3,i3,a,2e17.8)') 'ak0(', i, j, ') = ', ak0(i, j)
                     write (*, '(a,2e17.8)') 'ak1-ak0', ak1(i, j) - ak0(i, j)
                     print *, 'calc error'
                     pause
                     stop
                  end if
               end do
            end do
         else
            do i = 1, nz
               do j = 1, nk
                  if (ak1(i, j) /= ak0(i, j)) then
                     write (*, '(/)')
                     print *, 'z = ', i
                     print *, 'x = ', j
                     print *, 't = ', it
                     write (*, '(a,i3,i3,a,2e17.8)') 'a1(', i, j, ') = ', ak1(i, j)
                     write (*, '(a,i3,i3,a,2e17.8)') 'a0(', i, j, ') = ', ak0(i, j)
                     write (*, '(a,2e17.8)') 'ak1-ak0', ak1(i, j) - ak0(i, j)
                     print *, 'calc error'
                     pause
                     stop
                  end if
               end do
            end do
         end if
         !print *, 'perenos - check - ok'
         !for the next step
         jk0 = jk1
         ak0 = ak1
         a0 = a1
             !!dir$ else
         !print *, 'ok'
      else
         !for the next step
         jk0 = jk1
         ak0 = ak1
         a0 = a1
            !!dir$ endif
      end if

      !axt_z0(:,it+1) = a1(1,:)
      !axt_lz(:,it+1) = a1(nz,:)
      istat = fseek(222, 0, seek_end)
      inquire (222, nextrec=nr) !move the pointer to the end of the file (beginning of the last record)
      nr = nr + 1 !move the pointer to the end of the file (for the last record)

      write (222, rec=nr) ak1(nz, :)
      write (8, '(e17.8)', err=103) (it - 1)*h

      !call fn_for_fortran_to_call(ptr) ! calling a piece from c++

      if (recount == .false.) then
         write (*, '(a,f17.8,a,f5.3,a,e17.8, a,e17.8, a,e17.8, a,\,a)') 'time = ', it*h, &
            '   a(', ix_out*hx, ') = ', cdabs(atmp(ix_out)), '   eff = ', sum_eff, '   omega = ', omega, '       c', char(13)
      else
         write (*, '(a,f17.8,a,f5.3,a,e17.8, a,e17.8, a,e17.8, a,\,a)') 'time = ', it*h, &
            '   a(', ix_out*hx, ') = ', cdabs(atmp(ix_out)), '   eff = ', sum_eff, '   omega = ', omega, '       r', char(13)
      end if
   end do do_t

   write (*, '(/,/)')

   !closing the files to record the results
   close (222)
   close (8)
   close (3)
   close (53)
   close (553)

   !for drawing on charts
   !a_amp_z0_tend = cdabs(a1_z0(1,:))
   !a_amp_zl_tend = cdabs(a1_zl(nz,:))
   !ak_amp_z0_tend = cdabs(fs(a1_z0(1,:)))
   !ak_amp_zl_tend = cdabs(fs(a1_zl(nz,:)))

   call cpu_time(finish_time)
   print *, 'computation time = ', finish_time - start_time, ' seconds'

   call write_result()
   call finish()
   !ires = delfilesqq('$$$.bin')

   print *, 'finish!'
   pause
   stop
101 print *, 'error of file open.'
   pause
   stop
102 print *, 'error of file reading.'
   pause
   stop
103 print *, 'error of file writing.'
   pause
   stop

contains

   function jf(th)
      implicit none

      real(c_double), intent(in) :: th(:)
      complex(c_double_complex) jf

      jf = 2.0d0/dble(nth)*sum(cdexp(-im1*th))
   end function jf

   function rhs(ak, th)

      use ic, only: atmp, ixe1, ixe2, nt

      implicit none
      complex(c_double_complex), intent(in) :: ak(:)
      real(c_double), intent(in) :: th(:, :)
      real(c_double), dimension(size(th, 1), size(th, 2)) :: rhs, rhs1
      integer ix
      integer, save :: n = 0
      real(c_double) :: m = 0, m1 = 0

      !rhs1(:,1) = dreal(sum(fk1 * ak) * cdexp(im1 * th(:,1)))
      !rhs1(:,2) = dreal(sum(fk2 * ak) * cdexp(im1 * th(:,2)))

      rhs(:, 1) = dreal(mysum(fk1(1:nk)*ak)*cdexp(im1*th(:, 1)))
      rhs(:, 2) = dreal(mysum(fk2(1:nk)*ak)*cdexp(im1*th(:, 2)))

      !atmp = ifs(ak)
      !
      !rhs(:,1) = dreal(atmp(ixe1) * cdexp(im1 * th(:,1)))
      !rhs(:,2) = dreal(atmp(ixe2) * cdexp(im1 * th(:,2)))

      !if (dabs(maxval(rhs1)) > m1) m1 = dabs(maxval(rhs1))
      !if (dabs(maxval(rhs)) > m) m = dabs(maxval(rhs))

      !if (n == 1) then
      !    open(113, file = 'test.dat')
      !endif
      !write(113, '(i,8e17.8)') n, atmp(ixe1), mysum(fk1 * ak), atmp(ixe2), mysum(fk2 * ak)
      !if (n == nt-100) then
      !    close(113)
      !endif

      !n = n + 1

   end function rhs

   function mysum(a)
      implicit none
      complex(c_double_complex), intent(in) :: a(:)
      complex(c_double_complex) :: mysum
      integer(c_int) i, n

      n = size(a)
      mysum = dcmplx(0)

      do i = n, 1, -1
         mysum = mysum + a(i)
      end do
   end function mysum

   function dmysum(a)
      implicit none
      real(c_double), intent(in) :: a(:)
      real(c_double) :: dmysum
      integer(c_int) i, n

      n = size(a)
      !dmysum = dcmplx(0)
      dmysum = 0.0d0

      do i = n, 1, -1
         dmysum = dmysum + a(i)
      end do
   end function dmysum
   !end subroutine calculate_fortran
end program elektron2d

subroutine init() bind(c, name='init')

   use, intrinsic :: iso_c_binding
   use ic, only: nz, nx, nt, th, dthdz, a0, a1, g, k2, start_time, calc_idx, calc_theta, hx, h, &
                 fk1, fk2, xe, recount, ak0, lz, lt, delta, omega_d, lx, cold, nk, atmp, ops
   use fourier

   implicit none

   integer i, j

   interface
      subroutine read_param() bind(c, name='read_param')
      end subroutine read_param
      subroutine write_param() bind(c, name='write_param')
      end subroutine write_param
      function f_fn() result(f_res)
         use ic
         real(c_double), dimension(nx) :: f_res
      end function f_fn
      function g_fn() result(asm_res)
         use ic
         real(c_double), dimension(nx) :: asm_res
      end function g_fn
      function k2_fn() result(k2_res)
         use ic
         complex(c_double_complex), dimension(nk) :: k2_res
      end function k2_fn
      function dn_fn() result(dn_res)
         use ic
         complex(c_double_complex), dimension(nk) :: dn_res
      end function dn_fn
      function fk_fn(xe) result(fk_res)
         use, intrinsic :: iso_c_binding, only: c_double
         use ic, only: nk
         real(c_double), dimension(nk) :: fk_res
         real(c_double) xe
      end function fk_fn
      subroutine at0_fn(a0)
         use, intrinsic :: iso_c_binding, only: c_double_complex
         use ic, only: nz, nx
         complex(c_double_complex), dimension(nz, nx) :: a0
      end subroutine at0_fn
   end interface

   call read_param()
   !call write_param()
   call calc_idx()
   call allocate_arrays()
   call calc_zxt()
   call sincost_init(nx)
   call fft_init(nz)

   delta = delta + omega_d !new delta

   !smooth f
   if (cold .eq. .false.) then
      fk1(:) = fk_fn(xe)
      fk2(:) = fk_fn(lx - xe)
   else
      fk1 = 0.0d0
      fk2 = 0.0d0
   end if

   !smoothing function for a
   g = g_fn()

   write (*, *) 'lz = ', lz
   write (*, *) 'nz = ', nz
   write (*, *) 'lt = ', lt
   write (*, *) 'nt = ', nt
   write (*, *) 'omega_d = ', omega_d
   write (*, *) 'new delta = ', delta

   !k**2
   if (ops .eq. .false.) then
      k2 = k2_fn()
   else
      k2 = dn_fn()
   end if

   call cpu_time(start_time)

   open (1, file='init.dat')
   do i = 1, nk
      write (1, '(2e17.8)') fk1(i), fk2(i)
   end do
   close (1)

   atmp = g
   call sint(atmp)
   open (1, file='initg.dat')
   do i = 1, nx
      write (1, '(4e17.8)') (i - 1)*hx, g(i), dreal(atmp(i)), cdabs(atmp(i))
   end do
   close (1)

   !call isint(atmp)
   !open(1,file='initg_test.dat')
   !do i=1,nx
   !    write(1,'(4e17.8)') (i-1)*hx, g(i), dreal(atmp(i)), dimag(atmp(i))
   !enddo
   !close(1)

   return
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'

end subroutine init

subroutine finish() bind(c, name='finish')
   use fourier, only: sincost_destroy, fft_destroy

   call sincost_destroy()
   call fft_destroy()
   call deallocate_arrays()
end subroutine finish

subroutine calc_zxt()
   use ic, only: nz, nx, nt, z, x, h, hx

   implicit none

   integer i

   do i = 1, nz
      z(i) = (i - 1)*h
   end do
   do i = 1, nx
      x(i) = (i - 1)*hx
   end do

   open (7, file='z.dat', err=101)
   do i = 1, nz
      write (7, '(e17.8)', err=103) z(i)
   end do
   close (7)

   open (6, file='x.dat', err=101)
   do i = 1, nx
      write (6, '(e17.8)', err=103) x(i)
   end do
   close (6)

   return
101 stop 'error of file open.'
103 stop 'error of file writing.'
end subroutine calc_zxt

function a0_fn(n_it) result(a0_res)
    use ic, only : h, nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, iimp_t0, in_type, iimp_tend, imp_t0, imp_tsp, lx, central_mirror, xcp, alfa
   use fourier

   implicit none

   complex(c_double_complex), dimension(nx) :: a0_res, c
   integer(c_int), intent(in) :: n_it
   integer i, j, ng, n, seed(2), icp, ix(nx)
   real, dimension(5) :: rc

   j = n_it + 1 !move

   !initial conditions = 0 for times less t0 and large t0+tsp
   if ((j < iimp_t0) .or. (j > iimp_tend)) then
      a0_res = 0.0d0
      return
   end if

   if (in_type == 1) then
      !initial conditions for a (one pulse in the middle)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      a0_res(nx/2 + 2:) = 0.0d0
      do i = iimp_x0, nx/2 + 1
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0)) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 2) then
      !initial conditions for a (symmetric pulses at the edges)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
      do i = iimp_x0, iimp_xend
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0)) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 3) then
      !initial conditions from harmonics
      c = 0
      !seed = (/2147483562, 2147483398/)
      seed = (/3, 2/)
      call random_seed(size=n)
      if (n /= 2) stop 'error of random at a0_fn_stat'
      call random_seed(put=seed)
      call random_number(rc)
      !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

      c(3) = a0_peak

      a0_res = ifs(c)
      a0_res = cmplx(dreal(a0_res), 0.0d0) &
               *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))

      !a0_res = cmplx(dreal(a0_res),0.0d0)
      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 4) then
      !test initial conditions for a (one pulse in the middle)
      do i = 1, nx
         a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
   elseif (in_type == 6) then
      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)

         a0_res = a0_res*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      else
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx - 1/)

         a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)

         a0_res = a0_res*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end if
   else
      print *, 'error: wrong in_type'
      pause
      stop
   end if
end function a0_fn

function ak0_fn(n_it) result(ak0_res)
        use ic, only : h, nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, iimp_t0, in_type, iimp_tend, imp_t0, imp_tsp, lx, central_mirror, xcp, alfa, nk
   use fourier

   implicit none

   complex(c_double_complex), dimension(nk) :: ak0_res
   complex(c_double_complex), dimension(nx) :: a0_res, c
   integer(c_int), intent(in) :: n_it
   integer i, j, ng, n, seed(2), icp, ix(nx)
   real, dimension(5) :: rc

   j = n_it + 1 !move

   !initial conditions = 0 for times less t0 and large t0+tsp
   if ((j < iimp_t0) .or. (j > iimp_tend)) then
      a0_res = 0.0d0
      return
   end if

   if (in_type == 1) then
      !initial conditions for a (one pulse in the middle)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      a0_res(nx/2 + 2:) = 0.0d0
      do i = iimp_x0, nx/2 + 1
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0)) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 2) then
      !initial conditions for a (symmetric pulses at the edges)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
      do i = iimp_x0, iimp_xend
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0)) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 3) then
      !initial conditions from harmonics
      c = 0
      !seed = (/2147483562, 2147483398/)
      seed = (/3, 2/)
      call random_seed(size=n)
      if (n /= 2) stop 'error of random at a0_fn_stat'
      call random_seed(put=seed)
      call random_number(rc)
      !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

      c(3) = a0_peak

      a0_res = ifs(c)
      a0_res = cmplx(dreal(a0_res), 0.0d0) &
               *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))

      !a0_res = cmplx(dreal(a0_res),0.0d0)
      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 4) then
      !test initial conditions for a (one pulse in the middle)
      do i = 1, nx
         a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx) &
                     *dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end do
   elseif (in_type == 6) then
      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)

         a0_res = a0_res*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      else
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx - 1/)

         a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)

         a0_res = a0_res*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))*dsin(pi/imp_tsp*((j - 1)*h - imp_t0))
      end if
   else
      print *, 'error: wrong in_type'
      pause
      stop
   end if

   call sint(a0_res)
   ak0_res = a0_res(1:nk)
end function ak0_fn

function a0_fn_stat(n_it) result(a0_res)
   use ic, only: nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, in_type, lx, central_mirror, xcp, alfa!, coeff
   use fourier

   implicit none

   complex(c_double_complex), dimension(nx) :: a0_res, c
   real(c_double), dimension(nx) :: a0env
   integer(c_int), intent(in) :: n_it
   integer i, ng, n, seed(2), icp, ix(nx)
   real, dimension(5) :: rc
   real rn

   if (in_type == 1) then
      !initial conditions for a (one pulse in the middle)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      a0_res(nx/2 + 2:) = 0.0d0
      do i = iimp_x0, nx/2 + 1
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 2) then
      !initial conditions for a (symmetric pulses at the edges)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
      do i = iimp_x0, iimp_xend
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 3) then
      !initial conditions from harmonics
      c = 0
      !seed = (/2147483562, 2147483398/)
      !seed = (/3, 2/)
      !call random_seed(size = n)
      !if (n /= 2) stop 'error of random at a0_fn_stat'
      !call random_seed (put = seed)
      !call random_number(rc)
      !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

      c(2) = 0.1
      c(3) = 0.1
      c(4) = 0.05
      c(5) = 0.05

      !print *, size(c(2:10:2))
      !do i=1,size(rc)
      !    write(*,'(a,i2,a,f6.4,a,i2,a,f6.4,a,f6.4)') 'rc(', i, ') = ', rc(i), '     c(', 2*i, ') = ', dreal(c(2*i)), '   ', dimag(c(2*i))
      !enddo

      a0_res = ifs(c)
      a0_res = cmplx(dreal(a0_res), 0.0d0)

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 4) then
      !test initial conditions for a (one pulse in the middle)
      do i = 1, nx
         a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx)
      end do
   elseif (in_type == 5) then
      !specialist. initial conditions for a

      c = dcmplx(0)

      c(2) = dcmplx(0.1)
      !c(4) = dcmplx(0.1)
      !c(6) = dcmplx(0.1)
      !c(8) = dcmplx(0.1)
      !c(10) = dcmplx(0.1)
      !c(12) = dcmplx(0.1)
      !c(14) = dcmplx(0.1)

      a0_res = ifs(c)

      !open(1, file = 'test.dat')
      !do i=1,nx
      ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 6) then
      !specialist. initial conditions for a
      !initial conditions for a (symmetric pulses at the edges)

      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)
      else
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx - 1/)

         a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
      end if

      !open(1, file = 'test.dat')
      !do i=1,nx
      ! write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      ! !write(1, *) ix(i)
      !enddo
      !close(1)
      !stop
   elseif (in_type == 7) then
      !specialist. initial conditions for a
      !initial conditions for a (symmetric pulses at the edges)

      c = dcmplx(0)

      !c(2) = dcmplx(0.1)
      c(4) = dcmplx(0.1)
      c(6) = dcmplx(0.1)
      c(8) = dcmplx(0.1)
      c(10) = dcmplx(0.1)
      !c(12) = dcmplx(0.1)
      !c(14) = dcmplx(0.1)

      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0env(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0env(nx/2 + 2:nx) = a0env(nx/2:2:-1)
      else
         ix = (/1:nx/) - 1

         a0env = dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
      end if

      a0_res = ifs(c)*a0env

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !    !write(1, *) ix(i)
      !enddo
      !close(1)
      !stop
   else
      print *, 'error: wrong in_type'
      pause
      stop
   end if
end function a0_fn_stat

function ak0_fn_stat(n_it) result(ak0_res)
   use ic, only: nx, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, in_type, lx, central_mirror, xcp, alfa, nk
   use fourier

   implicit none

   complex(c_double_complex), dimension(nk) :: ak0_res
   complex(c_double_complex), dimension(nx) :: a0_res, c
   real(c_double), dimension(nx) :: a0env
   integer(c_int), intent(in) :: n_it
   integer i, ng, n, seed(2), icp, ix(nx)
   real, dimension(5) :: rc
   real rn

   if (in_type == 1) then
      !initial conditions for a (one pulse in the middle)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      a0_res(nx/2 + 2:) = 0.0d0
      do i = iimp_x0, nx/2 + 1
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 2) then
      !initial conditions for a (symmetric pulses at the edges)
      if (iimp_x0 > 1) a0_res(1:iimp_x0 - 1) = 0.0d0
      if (iimp_xend < nx) a0_res(iimp_xend + 1:nx) = 0.0d0
      do i = iimp_x0, iimp_xend
         a0_res(i) = a0_peak*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))*dsin(pi/imp_xsp*((i - 1)*hx - imp_x0))
      end do
      a0_res(nx:nx/2 + 2:-1) = a0_res(2:nx/2) !reflect pulse
   elseif (in_type == 3) then
      !initial conditions from harmonics
      c = 0
      !seed = (/2147483562, 2147483398/)
      !seed = (/3, 2/)
      !call random_seed(size = n)
      !if (n /= 2) stop 'error of random at a0_fn_stat'
      !call random_seed (put = seed)
      !call random_number(rc)
      !c(2:10:2) = cmplx(0.1 * rc(1:5), 0.0d0)

      c(2) = 0.1
      c(3) = 0.1
      c(4) = 0.05
      c(5) = 0.05

      !print *, size(c(2:10:2))
      !do i=1,size(rc)
      !    write(*,'(a,i2,a,f6.4,a,i2,a,f6.4,a,f6.4)') 'rc(', i, ') = ', rc(i), '     c(', 2*i, ') = ', dreal(c(2*i)), '   ', dimag(c(2*i))
      !enddo

      a0_res = ifs(c)
      a0_res = cmplx(dreal(a0_res), 0.0d0)

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i)), abs(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 4) then
      !test initial conditions for a (one pulse in the middle)
      do i = 1, nx
         a0_res(i) = a0_peak*dsin(1*pi/lx*(i - 1)*hx)
      end do
   elseif (in_type == 5) then
      !specialist. initial conditions for a

      c = dcmplx(0)

      c(2) = dcmplx(0.1)
      !c(4) = dcmplx(0.1)
      !c(6) = dcmplx(0.1)
      !c(8) = dcmplx(0.1)
      !c(10) = dcmplx(0.1)
      !c(12) = dcmplx(0.1)
      !c(14) = dcmplx(0.1)

      a0_res = ifs(c)

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !enddo
      !close(1)
      !stop
   elseif (in_type == 6) then
      !specialist. initial conditions for a
      !initial conditions for a (symmetric pulses at the edges)

      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0_res(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0_res(nx/2 + 2:nx) = a0_res(nx/2:2:-1)
      else
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx - 1/)

         a0_res = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
      end if

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !    !write(1, *) ix(i)
      !enddo
      !close(1)
      !stop
   elseif (in_type == 7) then
      !specialist. initial conditions for a
      !initial conditions for a (symmetric pulses at the edges)

      c = dcmplx(0)

      !c(2) = dcmplx(0.1)
      c(4) = dcmplx(0.1)
      c(6) = dcmplx(0.1)
      c(8) = dcmplx(0.1)
      c(10) = dcmplx(0.1)
      !c(12) = dcmplx(0.1)
      !c(14) = dcmplx(0.1)

      if (central_mirror == .false.) then
         icp = xcp/hx + 1
         iimp_xend = 2*icp - 1
         ix = 0; ix = (/0:nx/2 - 1/)

         a0env(1:nx/2) = a0_peak*dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
         a0env(nx/2 + 2:nx) = a0env(nx/2:2:-1)
      else
         ix = (/1:nx/) - 1

         a0env = dexp(-(ix*hx - xcp)**2/alfa) !+ dexp(-(ix * hx - xcp**2)**2/alfa)
      end if

      a0_res = ifs(c)*a0env

      !open(1, file = 'test.dat')
      !do i=1,nx
      !    write(1, '(4e17.8)') (i-1)*hx, dreal(a0_res(i)), dimag(a0_res(i))
      !    !write(1, *) ix(i)
      !enddo
      !close(1)
      !stop
   else
      print *, 'error: wrong in_type'
      pause
      stop
   end if

   call sint(a0_res)
   ak0_res = a0_res(1:nk)
end function ak0_fn_stat

subroutine at0_fn(a0)
   use ic, only: lz, h, nx, nz, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, iimp_t0, in_type, &
                 iimp_tend, imp_t0, imp_tsp, lx, central_mirror, xcp, alfa
   use fourier

   implicit none

   interface
      function a0_fn(it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nx) :: a
         integer(c_int), intent(in) :: it
      end function a0_fn
   end interface

   complex(c_double_complex), dimension(nz, nx) :: a0
   integer i, j

   do j = 1, nz
      do i = 1, nx
         !a0(j,i) = 0.01 * dsin(1.0 * pi / lx * (i-1)*hx) * dsin(1.0 * pi / lx * (i-1)*hx) &
         !    * dsin(1.0 * pi / lz * (j-1)*h) * dsin(1.0 * pi / lz * (j-1)*h)
         !a0(j,i) = a0_peak * dsin(2.0 * pi / lx * (i-1)*hx) * dsin(2.0 * pi / lx * (i-1)*hx) &
         !    * dsin(1.0 * pi / lz * (j-1)*h) * dsin(1.0 * pi / lz * (j-1)*h)
         a0(j, :) = 0.0d0
      end do
      !a0(j,:) = a0_fn(j)
   end do
end subroutine at0_fn

subroutine akt0_fn(ak0)
   use ic, only: lz, h, nx, nz, hx, iimp_x0, a0_peak, pi, imp_xsp, imp_x0, iimp_xend, iimp_t0, in_type, &
                 iimp_tend, imp_t0, imp_tsp, lx, central_mirror, xcp, alfa, nk
   use fourier

   implicit none

   interface
      function a0_fn(it) result(a)
         use, intrinsic :: iso_c_binding
         use ic
         complex(c_double_complex), dimension(nx) :: a
         integer(c_int), intent(in) :: it
      end function a0_fn
   end interface

   complex(c_double_complex), dimension(nz, nk) :: ak0
   integer i, j

   !do j=1,nz
   !    do i=1,nx
   !        !a0(j,i) = 0.01 * dsin(1.0 * pi / lx * (i-1)*hx) * dsin(1.0 * pi / lx * (i-1)*hx) &
   !        !    * dsin(1.0 * pi / lz * (j-1)*h) * dsin(1.0 * pi / lz * (j-1)*h)
   !        !a0(j,i) = a0_peak * dsin(2.0 * pi / lx * (i-1)*hx) * dsin(2.0 * pi / lx * (i-1)*hx) &
   !        !    * dsin(1.0 * pi / lz * (j-1)*h) * dsin(1.0 * pi / lz * (j-1)*h)
   !        a0(j,:) = 0.0d0
   !    enddo
   !    !a0(j,:) = a0_fn(j)
   !enddo

   ak0 = dcmplx(0.0d0)
end subroutine akt0_fn

function fk_fn(xe) result(fk_res)

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use ic, only: nk, pi, hx, lx

   implicit none

   real(c_double) :: fk_res(nk), xe!, arg(nx)
   integer(c_int) n(nk), i

   n = (/0:nk - 1/)

   !arg = pi * n  * xe / lx
   fk_res = dsin(pi*n*xe/lx)

   !if (xe < lx/2) then
   !    open(1, file = 'test1.dat')
   !    do i=1,nx
   !        write(1,'(i,2e17.8)') i-1, fk_res(i), arg(i)
   !    enddo
   !    close(1)
   !else
   !    open(1, file = 'test2.dat')
   !    do i=1,nx
   !        write(1,'(i,6e17.8,i)') i-1, fk_res(i), arg(i), xe, lx, xe/lx, pi, n(i)
   !    enddo
   !    close(1)
   !endif

end function fk_fn

!function a0_fn_stat() result(a0_res)
! use ic
!
! implicit none
!
! complex(c_double_complex), dimension(nx) :: a0_res
! integer i
!
! !initial conditions for a (one pulse in the middle)
! if (iimp_x0 > 1) a0_res(1:iimp_x0-1) = 0.0d0
! a0_res(nx/2+2:) = 0.0d0
! do i=iimp_x0,nx/2+1
! a0_res(i) = a0_peak * dsin(pi / imp_xsp * ((i-1)*hx - imp_x0)) * dsin(pi / imp_xsp * ((i-1)*hx - imp_x0))
! end do
! a0_res(nx:nx/2+2:-1) = a0_res(2:nx/2) !reflect pulse
!
! !initial conditions for a (symmetric pulses at the edges)
! !if (iimp_x0 > 1) a0_res(1:iimp_x0-1) = 0.0d0
! !if (iimp_end < nx) a0_res(iimp_end+1:nx) = 0.0d0
! !do i=iimp_x0,iimp_end
! ! a0_res(i) = a0_peak * dsin(pi / imp_xsp * ((i-1)*hx - imp_x0)) * dsin(pi / imp_xsp * ((i-1)*hx - imp_x0))
! !end do
! !a0_res(nx:nx/2+2:-1) = a0_res(2:nx/2) !reflect pulse
!end function a0_fn_stat

!function f_fn() result(fsm)
! use ic
!
! implicit none
!
! real(c_double), dimension(nx) :: fsm
! integer i
!
! !smooth f
! if(ifsm_x0 > 1) fsm(1:ifsm_x0-1) = 0.0d0
! if (ifsm_end < nx/2+1) fsm(ifsm_end+1:nx/2+1) = fsm_amp
! if (ifsm_x0 /= ifsm_end) then
! do i=ifsm_x0,ifsm_end
! fsm(i) = fsm_amp * dsin(pi / 2.0d0 / fsm_sp * ((i-1)*hx - fsm_x0)) * dsin(pi / 2.0d0 / fsm_sp * ((i-1)*hx - fsm_x0))
! end do
! else
! fsm(ifsm_end) = fsm_amp
! end if
! fsm(nx:nx/2+2:-1) = fsm(2:nx/2) !reflect f
!end function f_fn

function g_fn() result(g_res)
   use ic

   implicit none

   real(c_double), dimension(nx) :: g_res
   integer i

   ig0 = g_x0/hx + 1
   ig1 = g_x1/hx + 2

   if (central_mirror == .false.) then
      g_res = 1.0d0
      g_res(ig0:ig1) = 0.0d0
   else
      g_res = 0.0d0
      g_res(ig0:ig1) = 1.0d0
   end if

   g_res = g_res*g_amp

   !print *, i0, nx-i0+2,  i1
   !print *, 'apsize_minus_mirr = ', apsize_minus_mirr, sum(g_res)
   !pause
   !stop
end function g_fn
!function g_fn() result(asm_res)
! use ic
!
! implicit none
!
! real(c_double), dimension(nx) :: asm_res
! integer i
!
! !smoothing function for a
! if (iasm_x0 > 1) asm_res(1:iasm_x0-1) = 1.0d0
! if (iasm_xend < nx) asm_res(iasm_xend+1:nx) = 0.0d0
! if (iasm_x0 /= iasm_xend) then
! do i=iasm_x0,iasm_xend
! asm_res(i) = g_amp * dsin(pi / 2.0d0 + pi / 2.0d0 / g_x1 * ((i-1)*hx - g_x0)) * dsin(pi / 2.0d0 + pi / 2.0d0 / g_x1 * (( i-1)*hx - g_x0))
! end do
! else
! asm_res(iasm_xend) = 0.0
! end if
! asm_res(nx:nx/2+2:-1) = asm_res(2:nx/2) !reflect
!end function g_fn

function k2_fn() result(k2_res)
   use ic

   implicit none

   complex(c_double_complex), dimension(nk) :: k2_res
   integer i
   real(c_double) w

   !k**2
   do i = 1, nk
      w = pi*(i - 1)/lx
      !k2_res(i) = w * w - was
      !k2_res(i) = - w * w ! become
      k2_res(i) = -w*w/kk! become
   end do

   open (1, file='k2_n.dat')
   do i = 1, nk
      write (1, '(i,2e17.8)') i, k2_res(i)
   end do
   close (1)
end function k2_fn

function dn_fn() result(dn_res)
   use ic, only: nk, c_double_complex, c_double, lambda, lx, im1, pi

   implicit none

   complex(c_double_complex), dimension(nk) :: dn_res
   complex(c_double_complex) k
   real(c_double) tmp
   integer i

   k = 2.0d0*pi/lambda

   dn_res(1) = dcmplx(1)
   do i = 1, nk
      tmp = 1.0d0 - (i - 1)*(i - 1)/4.0d0*(lambda/lx)*(lambda/lx)
      dn_res(i) = dsqrt(tmp) - 1.0d0
   end do

   dn_res = k*dn_res

   !tmp = 1.0d0 - 1.0d0 / 4.0d0 * (lambda / lx) * (lambda / lx)
   !if (tmp .ge. 0.0d0) then
   !    dn1 = dsqrt(tmp)
   !else
   !    dn1 = im1 * dsqrt(dabs(tmp))
   !endif
   !
   !do i=1,nk
   !    tmp = 1.0d0 - (i * i) / 4.0d0 * (lambda / lx) * (lambda / lx)
   !    if (tmp .ge. 0.0d0) then
   !        dn_res(i) = dsqrt(tmp) - dn1
   !    else
   !        dn_res(i) = im1 * dsqrt(dabs(tmp)) - dn1
   !    endif
   !end do

   open (1, file='delta_n.dat')
   do i = 1, nk
      write (1, '(i,2e17.8)') i, dn_res(i)
   end do
   close (1)
end function dn_fn

subroutine allocate_arrays()
   use ic

   implicit none

   integer(c_int) err_alloc

   allocate (a0(nz, nk), a1(nz, nk), ak0(nz, nk), ak1(nz, nk), jk1(nz, nk), jk0(nz, nk), atmp(nx), a0z(nx), aktmp(nx), &
             ak_amp_z0_tend(nx), ak_amp_zl_tend(nx), &
             th(nth, nz, 2), dthdz(nth, nz, 2), fk1(nk), fk2(nk), rhs0(nth, 2), z(nz), x(nx), k2(nk), &
             a_amp_z0_tend(nx), a_amp_zl_tend(nx), g(nx), ex(nk), dlt(nk), &
        akz0(nk,nz), jkz0(nk,nz), ak0z0(nk), a0z0(nx), sum_abs2_a_plus_by_z(nx), sum_abs2_a_plus_by_z_k(nk), a0z0cut(nx), azl(nx), akzl(nk), stat=err_alloc)   

   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if
end subroutine allocate_arrays

subroutine deallocate_arrays()
   use ic

   implicit none

   integer(c_int) err_dealloc

   deallocate (a0, a1, ak0, ak1, jk1, jk0, atmp, a0z, aktmp, &
               ak_amp_z0_tend, ak_amp_zl_tend, &
               th, dthdz, fk1, fk2, rhs0, z, x, k2, &
               a_amp_z0_tend, a_amp_zl_tend, g, ex, dlt, &
               akz0, jkz0, ak0z0, a0z0, sum_abs2_a_plus_by_z, sum_abs2_a_plus_by_z_k, a0z0cut, azl, akzl, stat=err_dealloc)

   if (err_dealloc /= 0) then
      print *, "deallocation error"
      pause
      stop
   end if
end subroutine deallocate_arrays

subroutine read_param() bind(c, name='read_param')
   use ic

   implicit none

   open (unit=1, file='input_fortran.in', status='old', err=101)
   read (unit=1, nml=param, err=102)
   close (unit=1)

   write (*, nml=param)

   return
101 print *, "error of file open"; pause; stop
102 print *, 'error of reading file "input_fortran.in"'; pause; stop
end subroutine read_param

subroutine write_param() bind(c, name='write_param')
   use ic

   implicit none

   open (unit=1, file='input_fortran.in', err=101)
   !write(unit=1,nml=param, err=103)
   close (unit=1)

   write (*, nml=param)

   return
101 stop "error of file open"
103 stop 'error of file writing.'
end subroutine write_param

subroutine write_result()
   use ic
   use, intrinsic :: iso_c_binding

   implicit none

   integer(c_int) i

   call cpu_time(start_time)

   open (1, file='aend.dat', err=101)
   do i = 1, nx
      write (1, '(1p3e17.8)', err=103) (i - 1)*hx, a_amp_z0_tend(i), a_amp_zl_tend(i)
   end do
   close (1)

   call cpu_time(finish_time)
   print *, 'writing time = ', finish_time - start_time, ' seconds'

   return
101 stop 'error of file open.'
102 stop 'error of file reading.'
103 stop 'error of file writing.'
end subroutine write_result

!subroutine write_result(dir, length) bind(c,name='write_result')
! use ic
! use, intrinsic :: iso_c_binding
!
! implicit none
!
! integer(c_int), intent(in), value :: length
! character(c_char), dimension(*), intent(in) :: dir
! character(len=length) file_name
! integer(c_int) i
!
! call cpu_time(start_time)
!
! file_name = transfer(dir(1:length),file_name)
!
! open(3,file='f.dat')
! do i=1,nx
! write(3,'(2e17.8)') x(i), f(i)
! end do
! close(3)
!
! open(unit=1,file=file_name // '/input_fortran.in', err=101)
! write(unit=1,nml=param, err=102)
! close(unit=1)
!
! open(21, file=file_name // '/axtz0.bin', err=101, form='binary')
! tmpxt = cdabs(axt_z0)
! write(21,err=103) nx,nt,x,t,tmpxt
! close(21)
!
! open(15, file=file_name // '/axtzlbin', err=101, form='binary')
! tmpxt = cdabs(axt_lz)
! write(15,err=103) nx,nt,x,t,tmpxt
! close(15)
!
! open(6, file=file_name // '/x.dat', err=101)
! do i=1,nx
! write(6,'(e17.8)') x(i)
! end do
! close(6)
!
! open(7, file=file_name // '/z.dat', err=101)
! do i=1,nz
! write(7,'(e17.8)') z(i)
! end do
! close(7)
!
! open(3, file=file_name // '/out_ak.dat', err=101)
! do i=1,nz
! sum_abs2_a_at_zl(i) = sum(cdabs(a1(i,:))*cdabs(a1(i,:)))
! sumk(i) = sum(eff(i,:) * f)
! write(3, 104, err=103) (i-1)*h, sum_a_at_zl(i), sumk(i), sum_a_at_zl(i) - sum_a_at_zl(1) - 4.0 * sumk(i)
! end do
! close(3)
! 104 format(4e17.8)
!
! open(13, file=file_name // '/am.bin', err=101, form='binary')
! tmpzx=cdabs(a1)
! write(13,err=103) nz,nx,z,x,tmpzx
! close(13)
!
! open(5, file=file_name // '/km.bin', err=101, form='binary')
! write(5,err=103) nz,nx,z,x,eff
! close(5)
!
! open(11, file=file_name // '/amr.bin', err=101, form='binary')
! tmpzx = real(a1)
! write(11,err=103) nz,nx,z,x,tmpzx
! close(11)
!
! open(12, file=file_name // '/ami.bin', err=101, form='binary')
! tmpzx = dimag(a1)
! write(12,err=103) nz,nx,z,x,tmpzx
! close(12)
!
! call sint(a1(1,:))
!
! !to spectrum output
! do i=1,nz
! call sint(a1(i,:))
! end do
!
! open(8, file=file_name // '/am_spec.bin', err=101, form='binary')
! tmpzx=cdabs(a1)
! write(8,err=103) nz,nx,z,x,tmpzx
! close(8)
!
! open(9, file=file_name // '/am_spec_r.bin', err=101, form='binary')
! tmpzx = dreal(a1)
! write(9,err=103) nz,nx,z,x,tmpzx
! close(9)
!
! open(10, file=file_name // '/am_spec_i.bin', err=101, form='binary')
! tmpzx = dimag(a1)
! write(10,err=103) nz,nx,z,x,tmpzx
! close(10)
!
! call cpu_time(finish_time)
! print *, 'writing time = ', finish_time - start_time, 'seconds'
!
! call fftw_destroy()
!
! return
! 101 stop 'error of file open.'
! 102 stop 'error of file reading.'
! 103 stop 'error of file writing.'
!end subroutine write_result

!subroutine get_param(out_h, out_lz, out_lx, out_lt, out_nx, out_nth, out_delta, out_az0, out_imp_x0, out_imp_xsp, out_imp_t0, out_imp_tsp, &
! out_fsm_amp, out_fsm_x0, out_fsm_sp, out_r, out_asm_amp, out_asm_x0, out_asm_sp) bind(c,name='get_param')
! use ic
!
! implicit none
!
! real(c_double), intent(out) :: out_h, out_lz, out_lx, out_lt, out_delta, out_az0, out_imp_x0, out_imp_xsp, out_imp_t0, out_imp_tsp, &
! out_fsm_amp, out_fsm_x0, out_fsm_sp, out_r, out_asm_amp, out_asm_x0, out_asm_sp
! integer(c_int), intent(out) :: out_nx, out_nth
!
! out_h=h; out_lz=lz; out_lx = lx; out_lt=lt; out_nx = nx; out_nth = nth; out_delta = delta
! out_az0 = a0_peak; out_imp_x0 = imp_x0; out_imp_xsp = imp_xsp; out_imp_t0 = imp_t0; out_imp_tsp = imp_tsp; out_fsm_amp = fsm_amp;
! out_fsm_x0 = fsm_x0; out_fsm_sp = fsm_sp; out_r=r; out_asm_amp = g_amp; out_asm_x0 = g_x0; out_asm_sp = g_x1
!end subroutine get_param
!
!subroutine set_param(in_h, in_lz, in_lx, in_lt, in_nx, in_nth, in_delta, in_az0, in_imp_x0, in_imp_xsp, in_imp_t0, in_imp_tsp, in_fsm_amp, &
! in_fsm_x0, in_fsm_sp, in_r, in_asm_amp, in_asm_x0, in_asm_sp) bind(c,name='set_param')
! use ic
!
! implicit none
!
! real(c_double), intent(in) :: in_h, in_lz, in_lx, in_lt, in_delta, in_az0, in_imp_x0, in_imp_xsp, in_imp_t0, in_imp_tsp, in_fsm_amp, &
! in_fsm_x0, in_fsm_sp, in_r, in_asm_amp, in_asm_x0, in_asm_sp
! integer(c_int), intent(in) :: in_nx, in_nth
!
! h = in_h; lz=in_lz; lx = in_lx; lt = in_lt; nx=in_nx; nth = in_nth; delta=in_delta
! a0_peak = in_az0; imp_x0 = in_imp_x0; imp_xsp = in_imp_xsp; imp_t0 = in_imp_t0; imp_xsp = in_imp_xsp; fsm_amp = in_fsm_amp; fsm_x0 = in_fsm_x0;
! fsm_sp = in_fsm_sp; r = in_r; g_amp = in_asm_amp; g_x0 = in_asm_x0; g_x1 = in_asm_sp
!end subroutine set_param
!
!subroutine get_parameters_for_charts(px, paxz0, paxlz, pasp_amp_z0, pasp_amp_lz, out_nx) bind(c,name='get_parameters_for_charts')
! use, intrinsic :: iso_c_binding
! use ic
! use arrays
!
! implicit none
!
! type(c_ptr), intent(out) :: px, paxz0, paxlz, pasp_amp_z0, pasp_amp_lz
! integer(c_int), intent(out) :: out_nx
!
! px = c_loc(x)
! paxz0 = c_loc(a_amp_z0_tend)
! paxlz = c_loc(a_amp_zl_tend)
! pasp_amp_z0 = c_loc(ak_amp_z0_tend)
! pasp_amp_lz = c_loc(ak_amp_zl_tend)
! out_nx = nx
!end subroutine get_parameters_for_charts
!
!function fmin(a, n) result(m) bind(c,name='fmin')
! use, intrinsic :: iso_c_binding
!
! implicit none
!
! integer(c_int), intent(in) :: n
! real(c_double), dimension(n), intent(in) :: a
! real(c_double) m
!
! m = minval(a)
!end function fmin
!
!function fmax(a, n) result(m) bind(c,name='fmax')
! use, intrinsic :: iso_c_binding
!
! implicit none
!
! integer(c_int), intent(in) :: n
! real(c_double), dimension(n), intent(in) :: a
! real(c_double) m
!
! m = maxval(a)
!end function fmax

subroutine makea(atmp, ak)
   use ic, only: nk, nx
   use fourier
   use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

   implicit none

   complex(c_double_complex), dimension(:), intent(inout) :: atmp
   complex(c_double_complex), dimension(:), intent(in) :: ak
   integer(c_int) n1, n2

   n1 = size(ak)
   n2 = size(atmp)

   if (n1 .ne. nk .or. n2 .ne. nx) then
      print *, 'error in "makea"'
      pause
      stop
   end if

   atmp = dcmplx(0.0d0)
   atmp(1:nk) = ak

   call isint(atmp)
end subroutine makea

subroutine makeak(ak, atmp)
   use ic, only: nk, nx, aktmp
   use fourier
   use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

   implicit none

   complex(c_double_complex), dimension(:), intent(inout) :: ak
   complex(c_double_complex), dimension(:), intent(in) :: atmp
   integer(c_int) n1, n2

   n1 = size(ak)
   n2 = size(atmp)

   if (n1 .ne. nk .or. n2 .ne. nx) then
      print *, 'error in "makeak"'
      pause
      stop
   end if

   aktmp = atmp

   call sint(aktmp)
   ak = aktmp(1:nk)
end subroutine makeak

subroutine make_a0z0_ak1_atmp(a0z0, az0cut, ak0, ak)
   use ic, only: nk, nx, r0, g
   use fourier
   use, intrinsic :: iso_c_binding, only: c_double_complex, c_int

   implicit none

   interface
      subroutine makea(atmp, ak)
         use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
         complex(c_double_complex), dimension(:), intent(inout) :: atmp
         complex(c_double_complex), dimension(:), intent(in) :: ak
      end subroutine makea
      subroutine makeak(ak, atmp)
         use, intrinsic :: iso_c_binding, only: c_double_complex, c_int
         complex(c_double_complex), dimension(:), intent(inout) :: ak
         complex(c_double_complex), dimension(:), intent(in) :: atmp
      end subroutine makeak
   end interface

   complex(c_double_complex), dimension(:), intent(inout) :: a0z0, ak0, az0cut
   complex(c_double_complex), dimension(:), intent(in) :: ak
   integer(c_int) n1, n2, n3, n4

   n1 = size(ak)
   n2 = size(az0cut)
   n3 = size(a0z0)
   n4 = size(ak0)

   if (n1 .ne. nk .or. n2 .ne. nx .or. n3 .ne. nx .or. n4 .ne. nk) then
      print *, 'error in "makea"'
      pause
      stop
   end if

   call makea(a0z0, ak) !before cutting

   az0cut = a0z0*r0*g !mirror reflection in z=0 and cut

   call sint(az0cut)

   ak0 = az0cut(1:nk)

   call makea(az0cut, ak0)
end subroutine make_a0z0_ak1_atmp
