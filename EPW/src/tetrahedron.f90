module tetrahedron
USE kinds, ONLY : DP
implicit none
type point
real(DP) :: w(10),c(10),d(10),e(10)
integer :: i
end type point

type tetra
type(point) :: p(4)
end type tetra

end module tetrahedron

subroutine make_kp_reg(nx,ny,nz,wkt,ntet, tet)
! generate a regular mesh from 0 to g_i, with eventual shift of origin only
! Tetrahedron are generated. How are Tetra near BZ boundaries treated correctly?
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer :: i,j,k,l,nx,ny,nz,nk,nt,np,pt(8), ntet
real(DP) :: q(3),ep(3),kpt(3,nx*ny*nz),wkt(nx*ny*nz)
type(tetra) :: tet(ntet)

ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
! shft = (-0.5d0)*(g1+g2+g3)

do i=1,nx*ny*nz
	wkt(i)=1d0/(nx*ny*nz)
enddo

np = 0
nt = 0
do i = 1,nx
	do j = 1,ny
		do k = 1,nz			
			do l = 1,6
				pt(1) = (i-1)*ny*nz + (j-1)*nz + (k-1) +1
				pt(2) = pt(1)+1					
				pt(3) = pt(1)+nz
				pt(4) = pt(1)+nz+1
				pt(5) = pt(1)+nz*ny
				pt(6) = pt(1)+nz*ny+1
				pt(7) = pt(1)+nz*ny+nz
				pt(8) = pt(1)+nz*ny+nz+1
				if(k .eq. nz) then
					pt(2)=pt(2)-nz
					pt(4)=pt(4)-nz
					pt(6)=pt(6)-nz
					pt(8)=pt(8)-nz
				endif
				if(j .eq. ny) then
					pt(3)=pt(3)-ny*nz
					pt(4)=pt(4)-ny*nz
					pt(7)=pt(7)-ny*nz
					pt(8)=pt(8)-ny*nz
				endif
				if(i .eq. nx) then
					pt(5)=pt(5)-nx*ny*nz
					pt(6)=pt(6)-nx*ny*nz
					pt(7)=pt(7)-nx*ny*nz
					pt(8)=pt(8)-nx*ny*nz
				endif

				nt = nt+1
				np = np+1

				if (l.eq.1) then   			! 1,2,3,6
					tet(nt)%p(1)%i = pt(1)
					tet(nt)%p(2)%i = pt(2)
					tet(nt)%p(3)%i = pt(3)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.2) then			! 2,3,4,6
					tet(nt)%p(1)%i = pt(2)
					tet(nt)%p(2)%i = pt(3)
					tet(nt)%p(3)%i = pt(4)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.3) then			! 1,3,5,6
					tet(nt)%p(1)%i = pt(1)
					tet(nt)%p(2)%i = pt(3)
					tet(nt)%p(3)%i = pt(5)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.4) then			! 3,4,6,8
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(4)
					tet(nt)%p(3)%i = pt(6)
					tet(nt)%p(4)%i = pt(8)
				elseif (l.eq.5) then			! 3,5,6,7
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(5)
					tet(nt)%p(3)%i = pt(6)
					tet(nt)%p(4)%i = pt(7)
				else           			! 3,6,7,8
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(6)
					tet(nt)%p(3)%i = pt(7)
					tet(nt)%p(4)%i = pt(8)
				endif
			enddo	
		enddo
	enddo
enddo

! write(ulog,*)'KP_REG: Number of regular kpoints generated is=',nk
! 
! if (mod(nx,2).eq.0 .and. mod(ny,2).eq.0 .and. mod(nz,2).eq.0 ) then
! 	do i=1,nk
! 		kpt(:,i) = kpt(:,i)+shft(:)
! 	enddo
! endif
! 
! 2  format(i7,2x,3(1x,f12.5),5x,f9.5)
! 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
! close(126)

end subroutine make_kp_reg


subroutine eigen_tet(ntet,eival,tet,n,nkp)
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer j,k,l,ntet, n,nkp,w,cnd
real(DP), intent(in) :: eival(n,nkp)
type(tetra), intent(inout) :: tet(ntet)

do j=1,ntet !! this many tetra
	do k=1,4   !! corners
		do l=1,n   !! number of bands
			w = tet(j)%p(k)%i  !! i is the index of k points in the original kmesh
			tet(j)%p(k)%w(l)=eival(l,w) !! eigenvalues (band energies of band l, kpoint "w"
		enddo
	enddo
enddo

end subroutine eigen_tet

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eigen_tet_fermi(ntet,eival,tet,n,nkp)
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer j,k,l,ntet, n,nkp,w,cnd
real(DP), intent(in) :: eival(n,nkp)
type(tetra), intent(inout) :: tet(ntet)

do j=1,ntet !! this many tetra
	do k=1,4   !! corners
		do l=1,n   !! number of bands
			w = tet(j)%p(k)%i  !! i is the index of k points in the original kmesh
			tet(j)%p(k)%w(l)=eival(l,w) !! eigenvalues (band energies of band l, kpoint "w"
		enddo
	enddo
enddo

end subroutine eigen_tet_fermi

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine weight_tet(nktot, ntet,nb,e,tet,wkt)
!! input:  number of bands, the energy e that the weights are going to be computed, tet: tetrahedron
!! the first 5 parameters are easy. Just need to get tet right.
!! Lambin 1984
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer i,j,k,l,m,nb,nktot, ntet,uns(4),zero,kkk,ii
real(DP) :: e1,e2,e3,e4,e21,e31,e41,e32,e42,e43,en(4),en1(4),a1,a2,a3,a4,b11,b12,b21,b22,b31,b32,b41,b42
real(DP) :: a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44
real(DP) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,d41,d42,d43
real(DP) :: c11,c12,c21,c22,c31,c32,c41,c42
real(DP) :: s1,s2,s3,s4,s5,s6
real(DP) :: e,wkt(nktot)
type(tetra), intent(inout) :: tet(ntet)

zero=0

do i=1,ntet    ! allocate weights for imaginary terms
	do j=1,nb
		do k=1,4
			en(k)=tet(i)%p(k)%w(j)
		enddo
		en1=en
		call ssort(en,uns,4)  ! 
		if (en(1).ne.en(2) .and. en(2).ne.en(3) .and. en(3).ne.en(4)) then
			do k=1,4
				if (en(k).eq.en1(1)) then
					uns(k) = 1
				elseif (en(k).eq.en1(2)) then
					uns(k) = 2
				elseif (en(k).eq.en1(3)) then
					uns(k) = 3
				elseif (en(k).eq.en1(4)) then
					uns(k) = 4
				endif
			enddo
		else
			do k=1,4
				if (en(k).eq.en1(1)) then
					uns(k) = 1
				elseif (en(k).eq.en1(2)) then
					uns(k) = 2
				elseif (en(k).eq.en1(3)) then
					uns(k) = 3
				elseif (en(k).eq.en1(4)) then
					uns(k) = 4
				endif
				if (k.ne.1) then
					do l=1,k-1
						if (uns(l).eq.uns(k) .and. uns(k).eq.1) then
							if (en(k).eq.en1(2)) then
								uns(k) = 2
							elseif (en(k).eq.en1(3)) then
								uns(k) = 3
							elseif (en(k).eq.en1(4)) then
								uns(k) = 4
							endif
						elseif (uns(l).eq.uns(k) .and. uns(k).eq.2) then
							if (en(k).eq.en1(3)) then
								uns(k) = 3
							elseif (en(k).eq.en1(4)) then
								uns(k) = 4
							endif
						elseif (uns(l).eq.uns(k) .and. uns(k).eq.3) then
							uns(k) = 4
						endif
					enddo
				endif
			enddo  !! k = 1,4
		endif
		e1 = en(4)
		e2 = en(3)
		e3 = en(2)
		e4 = en(1)
		e21 = e2-e1
		e31 = e3-e1
		e41 = e4-e1
		e32 = e3-e2
		e42 = e4-e2
		e43 = e4-e3
! 		if(e21 .eq. 0. .or. e31 .eq. 0. .or. e41 .eq. 0. .or. e32 .eq. 0. .or. e42 .eq. 0. .or. e43 .eq. 0.) then 
! 			write(*,*) "Equal e tetra corners: ", i
! 		endif
		if (e.lt.e1 .or. e.gt.e4) then
			tet(i)%p(uns(4))%c(j) = 0
			tet(i)%p(uns(3))%c(j) = 0
			tet(i)%p(uns(2))%c(j) = 0
			tet(i)%p(uns(1))%c(j) = 0
		elseif (e1.le.e .and. e.le.e2) then
			tet(i)%p(uns(4))%c(j) = ((e2-e)/e21+(e3-e)/e31+(e4-e)/e41)*(e-e1)**2/(e41*e31*e21)
			tet(i)%p(uns(3))%c(j) = (e-e1)**3/(e21**2*e31*e41)
			tet(i)%p(uns(2))%c(j) = (e-e1)**3/(e21*e31**2*e41)
			tet(i)%p(uns(1))%c(j) = (e-e1)**3/(e21*e31*e41**2)
		elseif (e2.le.e .and. e.le.e3) then
			c11 = (e3-e)/e31**2
			c12 = (e4-e)/e41**2
			c21 = (e3-e)/e32**2
			c22 = (e4-e)/e42**2
			c31 = (e-e2)/e32**2
			c32 = (e-e1)/e31**2
			c41 = (e-e2)/e42**2
			c42 = (e-e1)/e41**2
			b11 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
			b12 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
			b21 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
			b22 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
			b31 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
			b32 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
			b41 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
			b42 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
			tet(i)%p(uns(4))%c(j) = .5*(c11*b11+c12*b12)
			tet(i)%p(uns(3))%c(j) = .5*(c21*b21+c22*b22)
			tet(i)%p(uns(2))%c(j) = .5*(c31*b31+c32*b32)
			tet(i)%p(uns(1))%c(j) = .5*(c41*b41+c42*b42)
		elseif (e3.le.e .and. e.le.e4) then
			tet(i)%p(uns(4))%c(j) = (e4-e)**3/(e41**2*e42*e43)
			tet(i)%p(uns(3))%c(j) = (e4-e)**3/(e41*e42**2*e43)
			tet(i)%p(uns(2))%c(j) = (e4-e)**3/(e41*e42*e43**2)
			tet(i)%p(uns(1))%c(j) = ((e-e3)/e43+(e-e2)/e42+(e-e1)/e41)*(e4-e)**2/(e41*e42*e43)
		endif
		tet(i)%p(uns(1))%c(j)=tet(i)%p(uns(1))%c(j)*wkt(tet(i)%p(uns(1))%i)*1.0d0/6.0d0 !! the weight is V_t/V_G which is the inverse of # tetra
		tet(i)%p(uns(2))%c(j)=tet(i)%p(uns(2))%c(j)*wkt(tet(i)%p(uns(2))%i)*1.0d0/6.0d0 !! spin=2 should not be included here, it's not relevant
		tet(i)%p(uns(3))%c(j)=tet(i)%p(uns(3))%c(j)*wkt(tet(i)%p(uns(3))%i)*1.0d0/6.0d0
		tet(i)%p(uns(4))%c(j)=tet(i)%p(uns(4))%c(j)*wkt(tet(i)%p(uns(4))%i)*1.0d0/6.0d0
		do ii=1,4
			if(isnan(tet(i)%p(uns(ii))%c(j) ) ) then
		 		tet(i)%p(uns(ii))%c(j)=0.0d0
		 	endif
		enddo
	enddo
enddo

end subroutine weight_tet

SUBROUTINE SSORT (X, IY, N)
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER N
REAL(DP) :: X(1:N)
INTEGER IY(N)
REAL(DP) :: TEMP
INTEGER I, ISWAP(1), ITEMP, ISWAP1
INTRINSIC MAXLOC
DO 200 I=1,N-1
	ISWAP=MAXLOC(X(I:N))
	ISWAP1=ISWAP(1)+I-1
	IF(ISWAP1.NE.I) THEN
		TEMP=X(I)
		X(I)=X(ISWAP1)
		X(ISWAP1)=TEMP
		ITEMP=IY(I)
		IY(I)=IY(ISWAP1)
		IY(ISWAP1)=ITEMP
	ENDIF
200 CONTINUE
RETURN
END
