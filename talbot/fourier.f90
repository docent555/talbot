module fourier
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'

   type(c_ptr), private :: plan_f, plan_b, plan_sc_f, plan_sc_b
   integer(c_int), private :: n_sc, n_fft, err_alloc
   complex(c_double_complex), private, allocatable :: in_f(:), out_f(:), in_b(:), out_b(:), &
                                                      in_sc_f(:), out_sc_f(:), in_sc_b(:), out_sc_b(:), &
                                                      w(:, :), tmp(:)
   complex(c_double_complex), parameter, private :: ic = (0.0d0, 1.0d0)
   real(c_double), private, parameter :: pi = 2.0d0*dacos(0.0d0)
   real(c_double), private :: l

contains
   subroutine dst_init(nn, ll)
      implicit none
      integer, intent(in) :: nn
      real(c_double), intent(in) :: ll
      real(c_double) dx
      integer i, j

      n_sc = nn
      l = ll

      dx = l/n_sc

      allocate (w(n_sc, n_sc), tmp(n_sc), stat=err_alloc)
      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      do i = 1, n_sc
         do j = 1, n_sc
            w(i, j) = dsin(pi*(i - 1)*(j - 1)*dx/l)
         end do
      end do
   end subroutine

   subroutine dst_destroy()
      deallocate (w)
   end subroutine dst_destroy

   subroutine dst(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)
      integer i, n

      n = n_sc/2
      tmp = dcmplx(0)

      do i = 1, n_sc
         tmp(i) = sum(v*w(i, :))/n
      end do

      v = tmp
   end subroutine dst

   subroutine idst(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)
      integer i

      tmp = dcmplx(0)

      do i = 1, n_sc
         tmp(i) = sum(v*w(:, i))
      end do

      v = tmp
   end subroutine idst

   subroutine fft_init(nn)
      implicit none
      integer, intent(in) :: nn

      n_fft = nn

      allocate (in_f(n_fft), out_f(n_fft), in_b(n_fft), out_b(n_fft), stat=err_alloc)
      if (err_alloc /= 0) stop "allocation error"

      plan_f = fftw_plan_dft_1d(n_fft, in_f, out_f, fftw_forward, fftw_estimate)
      plan_b = fftw_plan_dft_1d(n_fft, in_b, out_b, fftw_backward, fftw_estimate)
   end subroutine fft_init

   subroutine fft_destroy()
      call fftw_destroy_plan(plan_f)
      call fftw_destroy_plan(plan_b)
      deallocate (in_f)
      deallocate (out_f)
      deallocate (in_b)
      deallocate (out_b)
   end subroutine fft_destroy

   subroutine sincost_init(nn)
      implicit none
      integer, intent(in) :: nn

      n_sc = nn

      allocate (in_sc_f(2*n_sc), out_sc_f(2*n_sc), in_sc_b(2*n_sc), out_sc_b(2*n_sc), stat=err_alloc)
      if (err_alloc /= 0) stop "allocation error"

      plan_sc_f = fftw_plan_dft_1d(2*n_sc, in_sc_f, out_sc_f, fftw_forward, fftw_estimate)
      plan_sc_b = fftw_plan_dft_1d(2*n_sc, in_sc_b, out_sc_b, fftw_backward, fftw_estimate)
   end subroutine sincost_init

   subroutine sincost_destroy()
      call fftw_destroy_plan(plan_sc_f)
      call fftw_destroy_plan(plan_sc_b)
      deallocate (in_sc_f)
      deallocate (out_sc_f)
      deallocate (in_sc_b)
      deallocate (out_sc_b)
   end subroutine sincost_destroy

   subroutine fft(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)

      in_f = v

      call fftw_execute_dft(plan_f, in_f, out_f)

      v = out_f/n_fft
   end subroutine fft

   subroutine ifft(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)

      in_b = v

      call fftw_execute_dft(plan_b, in_b, out_b)

      v = out_b
   end subroutine ifft

   subroutine sint(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)

      in_sc_f(2:n_sc) = v(2:n_sc)
      in_sc_f(2*n_sc:n_sc + 2:-1) = -v(2:n_sc)

      in_sc_f(1) = (0.0d0, 0.0d0); in_sc_f(n_sc + 1) = (0.0d0, 0.0d0)

      call fftw_execute_dft(plan_sc_f, in_sc_f, out_sc_f)

      v = out_sc_f(1:n_sc)/(-ic*n_sc)
   end subroutine sint

   subroutine isint(v)
      implicit none
      complex(c_double_complex), intent(inout) :: v(:)

      in_sc_b(2:n_sc) = -v(2:n_sc)*(ic/2.0d0)
      in_sc_b(2*n_sc:n_sc + 2:-1) = -in_sc_b(2:n_sc)

      in_sc_b(1) = (0.0d0, 0.0d0); in_sc_b(n_sc + 1) = (0.0d0, 0.0d0)

      call fftw_execute_dft(plan_sc_b, in_sc_b, out_sc_b)

      v = out_sc_b(1:n_sc)
   end subroutine isint

   function fft_fun(v)
      implicit none
      complex(c_double_complex), intent(in) :: v(:)
      complex(c_double_complex), dimension(size(v)) :: fft_fun

      fft_fun = v
      call fft(fft_fun)
   end function fft_fun

   function ifft_fun(v)
      implicit none
      complex(c_double_complex), intent(in) :: v(:)
      complex(c_double_complex), dimension(size(v)) :: ifft_fun

      ifft_fun = v
      call ifft(ifft_fun)
   end function ifft_fun

   function fs(v)
      implicit none
      complex(c_double_complex), intent(in) :: v(:)
      complex(c_double_complex), dimension(size(v)) :: fs

      fs = v
      call sint(fs)

      !in_sc_f(2:n) = v(2:n)
      !in_sc_f(2*n:n+2:-1) = -v(2:n)
      !
      !in_sc_f(1) = (0.0d0,0.0d0); in_sc_f(n+1) = (0.0d0,0.0d0)
      !
      !call fftw_execute_dft(plan_sc_f, in_sc_f, out_sc_f)
      !
      !fs = out_sc_f(1:n)  / (-ic * n)
   end function fs

   function ifs(v)
      implicit none
      complex(c_double_complex), intent(in) :: v(:)
      complex(c_double_complex), dimension(size(v)) :: ifs

      ifs = v
      call isint(ifs)

      !in_sc_b(2:n) = -v(2:n) * (ic / 2.0d0)
      !in_sc_b(2*n:n+2:-1) = -in_sc_b(2:n)
      !
      !in_sc_b(1) = (0.0d0,0.0d0); in_sc_b(n+1) = (0.0d0,0.0d0)
      !
      !call fftw_execute_dft(plan_sc_b, in_sc_b, out_sc_b)
      !
      !ifs = out_sc_b(1:n)
   end function ifs
end module fourier
